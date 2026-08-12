"""Prior-predictive uncertainty propagation with resumable deterministic samples."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, Mapping, Optional

import numpy as np
import pandas as pd

from .analysis import _set_parameter
from .config import ModelConfig
from .model import simulate_0d


DEFAULT_PRIORS: Dict[str, Mapping[str, object]] = {
    "kinetics.capsule_release_s": {"distribution": "log_uniform", "low": 1e-7, "high": 3e-5},
    "kinetics.activity_multiplier": {"distribution": "uniform", "low": 0.5, "high": 5.0},
    "kinetics.k_oxygen_mol_m3": {"distribution": "log_uniform", "low": 0.005, "high": 0.20},
    "chemistry.calcite_rate_mol_m3_s": {"distribution": "log_uniform", "low": 1e-7, "high": 1e-2},
    "chemistry.wall_deposition_fraction": {"distribution": "uniform", "low": 0.05, "high": 1.0},
}
DEFAULT_SCENARIOS: Dict[str, Mapping[str, object]] = {
    "transport.crack_width_mm": {"distribution": "uniform", "low": 0.1, "high": 0.5},
    "environment.wet_hours_per_day": {"distribution": "uniform", "low": 6.0, "high": 24.0},
    "kinetics.agent_dosage_multiplier": {"distribution": "uniform", "low": 0.5, "high": 2.0},
}


def sample_prior(specification: Mapping[str, object], unit_value: float) -> float:
    """Transform a unit-interval value using a uniform or log-uniform prior."""

    low, high = float(specification["low"]), float(specification["high"])
    if high < low:
        raise ValueError("Prior upper bound must not be below lower bound")
    distribution = specification.get("distribution", "uniform")
    if distribution == "uniform":
        return low + float(unit_value) * (high - low)
    if distribution == "log_uniform":
        if low <= 0:
            raise ValueError("Log-uniform prior requires a positive lower bound")
        return float(np.exp(np.log(low) + float(unit_value) * (np.log(high) - np.log(low))))
    raise ValueError("Unsupported prior distribution: {}".format(distribution))


def _config_hash(config: ModelConfig) -> str:
    payload = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _outputs(config: ModelConfig, target_closure: float) -> Dict[str, float]:
    result = simulate_0d(config)
    frame = result.frame
    final = frame.iloc[-1]
    reached = frame.loc[frame["crack_closure_ratio"] >= target_closure, "time_d"]
    first_day = frame.loc[frame["time_d"] <= 1.0]
    initial_capsule = float(frame["capsule_calcium_lactate_mol_m3"].iloc[0])
    remaining_day1 = float(first_day["capsule_calcium_lactate_mol_m3"].iloc[-1])
    return {
        "closure_28d": float(final["crack_closure_ratio"]),
        "calcite_mass_mg": float(result.summary["calcite_mass_mg"]),
        "permeability_ratio": float(final["permeability_ratio"]),
        "transmissivity_ratio": float(final["transmissivity_ratio"]),
        "premature_consumption_day1": 1.0 - remaining_day1 / max(initial_capsule, 1e-30),
        "time_to_target_d": float(reached.iloc[0]) if len(reached) else np.nan,
    }


def prior_predictive(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    samples: int = 256,
    seed: int = 2026,
    priors: Mapping[str, Mapping[str, object]] = DEFAULT_PRIORS,
    scenarios: Mapping[str, Mapping[str, object]] = DEFAULT_SCENARIOS,
    target_closure: float = 0.5,
    resume: bool = False,
) -> Dict[str, object]:
    """Propagate epistemic priors and scenario variability without calibration claims."""

    base = copy.deepcopy(config or ModelConfig())
    base.simulation.days = 28.0
    base.simulation.output_interval_days = 1.0
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_path = output_dir / "prior_predictive_samples.csv"
    existing = pd.read_csv(sample_path) if resume and sample_path.exists() else pd.DataFrame()
    completed = set(existing["sample_id"].astype(int)) if len(existing) else set()
    rng = np.random.RandomState(seed)
    names = list(priors) + list(scenarios)
    unit_samples = rng.uniform(size=(max(int(samples), 0), len(names)))
    rows = existing.to_dict("records") if len(existing) else []
    for sample_id, unit_row in enumerate(unit_samples):
        if sample_id in completed:
            continue
        trial = copy.deepcopy(base)
        values = {}
        for index, name in enumerate(names):
            specification = priors[name] if name in priors else scenarios[name]
            value = sample_prior(specification, unit_row[index])
            _set_parameter(trial, name, value)
            values[name] = value
        rows.append({
            "sample_id": sample_id, **values, **_outputs(trial, target_closure),
            "evidence_class": "prior_predictive_model_output",
        })
        pd.DataFrame(rows).sort_values("sample_id").to_csv(sample_path, index=False)
    frame = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    metric_names = (
        "closure_28d", "calcite_mass_mg", "permeability_ratio",
        "transmissivity_ratio", "premature_consumption_day1", "time_to_target_d",
    )
    summaries = []
    for metric in metric_names:
        values = pd.to_numeric(frame.get(metric), errors="coerce").dropna()
        summaries.append({
            "metric": metric, "n": int(len(values)),
            "median": float(values.median()) if len(values) else np.nan,
            "lower_95": float(values.quantile(0.025)) if len(values) else np.nan,
            "upper_95": float(values.quantile(0.975)) if len(values) else np.nan,
            "failure_probability": float(np.mean(frame["closure_28d"] < target_closure)) if len(frame) else np.nan,
            "interval_type": "prior_predictive_interval",
        })
    pd.DataFrame(summaries).to_csv(output_dir / "prior_predictive_summary.csv", index=False)
    metadata = {
        "evidence_class": "prior_predictive_model_output",
        "interval_type": "prior_predictive_interval",
        "calibrated_prediction_interval": False,
        "epistemic_parameters": list(priors), "scenario_variables": list(scenarios),
        "samples": int(len(frame)), "random_seed": int(seed),
        "config_sha256": _config_hash(base), "resume_supported": True,
        "warning": "Intervals propagate literature/scenario priors and are not experimentally validated.",
    }
    (output_dir / "prior_predictive_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    return metadata

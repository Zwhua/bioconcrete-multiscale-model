"""Model-response counterfactual control and bottleneck analysis."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import hashlib
import json
from pathlib import Path
from typing import Callable, Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .analysis import _set_parameter
from .config import ModelConfig
from .evidence_state import UNCALIBRATED
from .model import simulate_0d


FACTORS: Mapping[str, str] = {
    "effective_activity": "kinetics.activity_multiplier",
    "substrate_inventory": "kinetics.capsule_calcium_lactate_mol_m3",
    "calcium_inventory": "chemistry.portlandite_mol_m3",
    "oxygen_transfer": "environment.oxygen_transfer_s",
    "capsule_release": "kinetics.capsule_release_s",
    "wall_deposition": "chemistry.wall_deposition_fraction",
    "crack_geometry": "transport.crack_width_mm",
    "csh_payload": "kinetics.capsule_csh_volume_fraction",
}
OUTPUTS = ("closure_28d", "calcite_mass_mg", "transmissivity_ratio")


def _get_parameter(config: ModelConfig, path: str) -> float:
    section, name = path.split(".", 1)
    return float(getattr(getattr(config, section), name))


def _outputs(config: ModelConfig) -> Dict[str, float]:
    result = simulate_0d(config)
    return {
        "closure_28d": float(result.summary["mean_crack_closure_ratio"]),
        "calcite_mass_mg": float(result.summary["calcite_mass_mg"]),
        "transmissivity_ratio": float(result.summary["mean_transmissivity_ratio"]),
    }


def control_coefficient(baseline: float, perturbed: float, relative_change: float) -> float:
    """Return a normalized local control coefficient, preserving undefined cases."""

    if relative_change == 0 or not np.isfinite(baseline) or abs(baseline) < 1e-30:
        return np.nan
    return float(((perturbed - baseline) / baseline) / relative_change)


def _task_id(factor: str, perturbation: float, prior_sample: int) -> str:
    value = "{}|{:.12g}|{}".format(factor, perturbation, prior_sample)
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]


def _config_from_dict(raw: Mapping[str, Mapping[str, object]]) -> ModelConfig:
    config = ModelConfig()
    for section, values in raw.items():
        for name, value in values.items():
            setattr(getattr(config, section), name, value)
    config.validate()
    return config


def _evaluate_task(payload: Tuple[Mapping[str, Mapping[str, object]], str, float, int, Mapping[str, float]]) -> Dict[str, object]:
    raw, factor, perturbation, prior_sample, baseline = payload
    base = _config_from_dict(raw)
    path = FACTORS[factor]
    value = _get_parameter(base, path)
    trial = copy.deepcopy(base)
    _set_parameter(trial, path, value * (1.0 + perturbation))
    trial.validate()
    changed = _outputs(trial)
    row: Dict[str, object] = {
        "task_id": _task_id(factor, perturbation, prior_sample),
        "prior_sample": prior_sample,
        "factor": factor,
        "parameter": path,
        "perturbation_fraction": perturbation,
        "baseline_parameter": value,
        "perturbed_parameter": value * (1.0 + perturbation),
        "evidence_label": UNCALIBRATED,
    }
    for output in OUTPUTS:
        row["baseline_{}".format(output)] = baseline[output]
        row["perturbed_{}".format(output)] = changed[output]
        row["control_{}".format(output)] = control_coefficient(
            baseline[output], changed[output], perturbation
        )
    return row


def summarize_controls(frame: pd.DataFrame, tie_fraction: float = 0.10) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Summarize nonlinear response and determine model-response bottlenecks."""

    rows = []
    for factor, selected in frame.groupby("factor"):
        record: Dict[str, object] = {"factor": factor}
        strengths = []
        for output in OUTPUTS:
            values = pd.to_numeric(selected["control_{}".format(output)], errors="coerce")
            negative = values[selected["perturbation_fraction"] < 0]
            positive = values[selected["perturbation_fraction"] > 0]
            neg = float(negative.median()) if negative.notna().any() else np.nan
            pos = float(positive.median()) if positive.notna().any() else np.nan
            record["negative_{}".format(output)] = neg
            record["positive_{}".format(output)] = pos
            record["median_abs_{}".format(output)] = float(values.abs().median())
            record["nonlinearity_{}".format(output)] = abs(pos - neg) if np.isfinite(pos + neg) else np.nan
            record["sign_reversal_{}".format(output)] = bool(np.isfinite(pos + neg) and pos * neg < 0)
            strengths.append(record["median_abs_{}".format(output)])
        record["aggregate_control"] = float(np.nanmax(strengths))
        rows.append(record)
    summary = pd.DataFrame(rows).sort_values("aggregate_control", ascending=False).reset_index(drop=True)
    if summary.empty or not np.isfinite(summary["aggregate_control"]).any():
        dominant, tied = "no_clear_dominant_bottleneck", []
    else:
        peak = float(summary["aggregate_control"].max())
        tied = summary.loc[summary["aggregate_control"] >= peak * (1.0 - tie_fraction), "factor"].tolist()
        dominant = tied[0] if len(tied) == 1 else ("parallel:" + ";".join(tied))
    metadata = {
        "dominant_bottleneck": dominant,
        "tied_factors": tied,
        "basis": "normalized counterfactual model responses",
        "evidence_label": UNCALIBRATED,
        "applicability_domain": "configured parameter priors and 0D 28-day model",
    }
    return summary, metadata


def counterfactual_bottleneck(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    perturbations: Sequence[float] = (-0.2, -0.1, 0.1, 0.2),
    workers: int = 1,
    resume: bool = False,
) -> Dict[str, object]:
    """Run deterministic, resumable model counterfactuals for each candidate control."""

    if not perturbations or any(value <= -1 or value == 0 for value in perturbations):
        raise ValueError("Perturbations must be nonzero and greater than -1")
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = output_dir / "counterfactual_control_coefficients.csv"
    existing = pd.read_csv(raw_path) if resume and raw_path.exists() else pd.DataFrame()
    rows = existing.to_dict("records") if len(existing) else []
    completed = set(existing.get("task_id", pd.Series(dtype=str)).astype(str))
    raw = (config or ModelConfig()).to_dict()
    baseline = _outputs(config or ModelConfig())
    tasks = [
        (raw, factor, float(delta), 0, baseline)
        for factor in FACTORS for delta in perturbations
        if _task_id(factor, float(delta), 0) not in completed
    ]
    failures = []

    def record(row: Dict[str, object]) -> None:
        rows.append(row)
        pd.DataFrame(rows).sort_values("task_id").to_csv(raw_path, index=False)

    if workers <= 1:
        for task in tasks:
            try:
                record(_evaluate_task(task))
            except Exception as error:
                failures.append({"task_id": _task_id(task[1], task[2], task[3]), "error": str(error)})
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_evaluate_task, task): task for task in tasks}
            for future in as_completed(futures):
                try:
                    record(future.result())
                except Exception as error:
                    task = futures[future]
                    failures.append({"task_id": _task_id(task[1], task[2], task[3]), "error": str(error)})
    frame = pd.DataFrame(rows).drop_duplicates("task_id", keep="last")
    frame.sort_values("task_id").to_csv(raw_path, index=False)
    controls, metadata = summarize_controls(frame)
    controls.to_csv(output_dir / "dominant_bottlenecks.csv", index=False)
    probability = controls[["factor", "aggregate_control"]].copy()
    total = probability["aggregate_control"].sum()
    probability["dominance_probability"] = probability["aggregate_control"] / total if total > 0 else np.nan
    probability["uncertainty_basis"] = "local perturbation envelope; prior ensemble pending"
    probability.to_csv(output_dir / "bottleneck_uncertainty.csv", index=False)
    metadata.update({"tasks": len(frame), "failures": failures, "workers": workers, "resume": resume})
    (output_dir / "bottleneck_summary.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata

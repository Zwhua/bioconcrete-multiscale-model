"""Mechanistic and baseline model-structure comparison interfaces."""

from __future__ import annotations

import copy
import json
import time
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .config import ModelConfig
from .model import simulate_0d


STRUCTURES = (
    "full_mechanistic", "no_biological_mineralization", "no_csh_fill",
    "no_environment_gate", "first_order_empirical", "zero_repair",
)


def information_criteria(observed: np.ndarray, predicted: np.ndarray, parameters: int) -> Dict[str, float]:
    """Return AIC/AICc only when real observations are explicitly supplied."""

    observed, predicted = np.asarray(observed, float), np.asarray(predicted, float)
    mask = np.isfinite(observed) & np.isfinite(predicted)
    n = int(mask.sum())
    if not n:
        return {"n": 0, "rmse": np.nan, "aic": np.nan, "aicc": np.nan}
    residual = predicted[mask] - observed[mask]
    rss = float(np.sum(residual**2))
    aic = float(n * np.log(max(rss / n, 1e-30)) + 2 * parameters)
    aicc = float(aic + 2 * parameters * (parameters + 1) / (n - parameters - 1)) if n > parameters + 1 else np.nan
    return {"n": n, "rmse": float(np.sqrt(rss / n)), "aic": aic, "aicc": aicc}


def compare_structures(
    output_dir: Path, config: Optional[ModelConfig] = None,
    observations: Optional[pd.DataFrame] = None,
) -> Dict[str, object]:
    """Compare model capabilities; omit error criteria when observations are absent."""

    base = copy.deepcopy(config or ModelConfig())
    rows = []
    parameter_counts = {
        "full_mechanistic": 6, "no_biological_mineralization": 4,
        "no_csh_fill": 5, "no_environment_gate": 4,
        "first_order_empirical": 2, "zero_repair": 0,
    }
    for structure in STRUCTURES:
        started = time.perf_counter()
        if structure in {"first_order_empirical", "zero_repair"}:
            closure = np.nan if structure == "first_order_empirical" else 0.0
            calcite, permeability = np.nan, 1.0
        else:
            trial = copy.deepcopy(base)
            if structure == "no_biological_mineralization":
                trial.kinetics.spore_density_rel = 0.0
                trial.kinetics.active_density_rel = 0.0
                trial.chemistry.calcite_rate_mol_m3_s = 0.0
            elif structure == "no_csh_fill":
                trial.kinetics.capsule_csh_volume_fraction = 0.0
            elif structure == "no_environment_gate":
                trial.kinetics.gate_logic = "static_suitability"
                trial.kinetics.activation_duration_h = 0.0
                trial.kinetics.response_delay_h = 0.0
            result = simulate_0d(trial)
            closure = result.summary["mean_crack_closure_ratio"]
            calcite = result.summary["calcite_mass_mg"]
            permeability = result.summary["mean_permeability_ratio"]
        metrics = {"n": 0, "rmse": np.nan, "aic": np.nan, "aicc": np.nan}
        if observations is not None and structure not in {"first_order_empirical"}:
            metrics = information_criteria(
                observations["crack_closure_ratio"].to_numpy(float),
                np.full(len(observations), closure), parameter_counts[structure],
            )
        rows.append({
            "structure": structure, "closure_28d": closure,
            "calcite_mass_mg": calcite, "permeability_ratio": permeability,
            "parameter_count": parameter_counts[structure], **metrics,
            "runtime_s": time.perf_counter() - started,
            "physically_interpretable": structure not in {"first_order_empirical", "zero_repair"},
            "answers_design_question": structure in {"full_mechanistic", "no_biological_mineralization", "no_csh_fill", "no_environment_gate"},
            "evidence_class": "model_structure_ablation",
        })
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "model_structure_comparison.csv", index=False)
    summary = {
        "structures": list(STRUCTURES), "observations_supplied": observations is not None,
        "information_criteria_available": observations is not None,
        "warning": "AIC/AICc are intentionally absent without observed data.",
    }
    (output_dir / "model_structure_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

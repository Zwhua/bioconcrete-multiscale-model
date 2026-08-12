"""D-optimal prospective experiment ranking for limited wet-lab capacity."""

from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd

from .config import ModelConfig


OBSERVABLE_PARAMETER_LINKS = {
    "closure_ratio": ("wall_deposition_fraction", "csh_payload", "calcite_rate"),
    "calcite_mass": ("calcite_rate", "calcium", "substrate"),
    "calcium_mM": ("calcium", "calcite_rate"),
    "substrate_mM": ("release_rate", "activity", "oxygen_limitation"),
}


def d_optimal_score(sensitivity: np.ndarray, prior_information: Optional[np.ndarray] = None) -> float:
    """Return log determinant gain for one candidate sensitivity matrix."""

    matrix = np.asarray(sensitivity, dtype=float)
    prior = np.eye(matrix.shape[1]) * 1e-6 if prior_information is None else np.asarray(prior_information)
    sign, value = np.linalg.slogdet(prior + matrix.T @ matrix)
    return float(value) if sign > 0 else -np.inf


def rank_experiments(output_dir: Path, config: Optional[ModelConfig] = None) -> Dict[str, object]:
    """Rank preregistered candidate measurements using a transparent D-opt proxy."""

    widths, doses, wettings = (0.1, 0.3, 0.5), (0.5, 1.0, 2.0), (6.0, 12.0, 24.0)
    times, observables = (1.0, 7.0, 14.0, 28.0), tuple(OBSERVABLE_PARAMETER_LINKS)
    parameter_names = sorted({p for values in OBSERVABLE_PARAMETER_LINKS.values() for p in values})
    rows = []
    for width, dose, wet, time_d, observable in itertools.product(widths, doses, wettings, times, observables):
        vector = np.zeros(len(parameter_names))
        for parameter in OBSERVABLE_PARAMETER_LINKS[observable]:
            vector[parameter_names.index(parameter)] = 1.0
        vector *= np.sqrt(time_d / 28.0) * (0.5 + wet / 24.0)
        if observable == "closure_ratio":
            vector *= 0.3 / width
        if observable in {"calcite_mass", "calcium_mM", "substrate_mM"}:
            vector *= dose
        score = d_optimal_score(vector.reshape(1, -1))
        experiment_id = hashlib.sha256(
            "{}|{}|{}|{}|{}".format(width, dose, wet, time_d, observable).encode()
        ).hexdigest()[:12]
        rows.append({
            "experiment_id": experiment_id, "crack_width_mm": width,
            "agent_dosage": dose, "wet_hours_per_day": wet, "time_d": time_d,
            "observable": observable, "d_optimal_score": score,
            "parameters_informed": ";".join(OBSERVABLE_PARAMETER_LINKS[observable]),
            "evidence_class": "model_informed_experimental_plan",
        })
    frame = pd.DataFrame(rows).sort_values("d_optimal_score", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    frame["plan"] = np.where(frame["rank"] <= 5, "minimum_executable",
                             np.where(frame["rank"] <= 10, "ideal_extension", "candidate"))
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "ranked_experiments.csv", index=False)
    frame.head(10).to_csv(output_dir / "recommended_experiments.csv", index=False)
    summary = {
        "method": "D-optimal local proxy", "candidate_count": len(frame),
        "minimum_executable_count": 5, "ideal_plan_count": 10,
        "status": "Model-informed experimental plan",
        "wet_lab_status": "Awaiting wet-lab execution",
        "experimental_validation": False,
    }
    (output_dir / "experiment_design_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

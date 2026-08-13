"""Time- and condition-resolved mechanistic model-structure comparison."""

from __future__ import annotations

import copy
import json
import time
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .config import ModelConfig
from .evidence import _predict_outputs
from .evidence_state import PUBLIC_CALIBRATED, UNCALIBRATED


STRUCTURES = (
    "full_mechanistic", "no_biological_mineralization", "no_csh_fill",
    "no_environment_gate", "first_order_empirical", "zero_repair",
)
PARAMETER_COUNTS = {
    "full_mechanistic": 6, "no_biological_mineralization": 4,
    "no_csh_fill": 5, "no_environment_gate": 4,
    "first_order_empirical": 2, "zero_repair": 0,
}


def information_criteria(observed: np.ndarray, predicted: np.ndarray, parameters: int) -> Dict[str, float]:
    """Return error and information criteria only for finite real observations."""

    observed, predicted = np.asarray(observed, float), np.asarray(predicted, float)
    mask = np.isfinite(observed) & np.isfinite(predicted)
    n = int(mask.sum())
    if not n:
        return {"n": 0, "mae": np.nan, "rmse": np.nan, "aic": np.nan, "aicc": np.nan,
                "prediction_interval_coverage": np.nan}
    residual = predicted[mask] - observed[mask]
    rss = float(np.sum(residual ** 2))
    aic = float(n * np.log(max(rss / n, 1e-30)) + 2 * parameters)
    aicc = float(aic + 2 * parameters * (parameters + 1) / (n - parameters - 1)) if n > parameters + 1 else np.nan
    return {
        "n": n, "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(rss / n)), "aic": aic, "aicc": aicc,
        "prediction_interval_coverage": np.nan,
    }


def _structure_config(base: ModelConfig, structure: str) -> ModelConfig:
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
    return trial


def _grouped_split(frame: pd.DataFrame) -> pd.Series:
    if "split" in frame and frame["split"].isin(["train", "internal_test"]).any():
        return frame["split"].fillna("train")
    specimens = frame.get("specimen_id", pd.Series(np.arange(len(frame)), index=frame.index)).astype(str)
    unique = sorted(specimens.unique())
    holdout = set(unique[::5]) if len(unique) >= 5 else set()
    return specimens.map(lambda value: "internal_test" if value in holdout else "train")


def _fit_empirical(train: pd.DataFrame) -> Optional[np.ndarray]:
    observed = pd.to_numeric(train["crack_closure_ratio"], errors="coerce").to_numpy(float)
    times = pd.to_numeric(train["time_d"], errors="coerce").to_numpy(float)
    mask = np.isfinite(observed) & np.isfinite(times)
    if mask.sum() < 4:
        return None
    fit = least_squares(
        lambda values: values[0] * (1.0 - np.exp(-values[1] * times[mask])) - observed[mask],
        np.array([min(max(np.nanmax(observed), 0.05), 1.0), 0.05]),
        bounds=([0.0, 1e-8], [1.0, 10.0]), max_nfev=200,
    )
    return fit.x


def compare_structures(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    observations: Optional[pd.DataFrame] = None,
    predictor: Callable[[ModelConfig, pd.DataFrame], pd.DataFrame] = _predict_outputs,
) -> Dict[str, object]:
    """Compare structures at every observation's actual time and condition."""

    base = copy.deepcopy(config or ModelConfig())
    output_dir.mkdir(parents=True, exist_ok=True)
    if observations is None:
        rows = [{
            "structure": structure, "parameter_count": PARAMETER_COUNTS[structure],
            "n": 0, "mae": np.nan, "rmse": np.nan, "aic": np.nan, "aicc": np.nan,
            "prediction_interval_coverage": np.nan, "runtime_s": 0.0,
            "answers_design_question": structure not in {"first_order_empirical", "zero_repair"},
            "evidence_label": UNCALIBRATED,
        } for structure in STRUCTURES]
        pd.DataFrame(rows).to_csv(output_dir / "model_structure_comparison.csv", index=False)
        summary = {"structures": list(STRUCTURES), "observations_supplied": False,
                   "information_criteria_available": False,
                   "warning": "AIC/AICc are intentionally absent without observed data."}
        (output_dir / "model_structure_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return summary

    frame = observations.reset_index(drop=True).copy()
    required = {"time_d", "crack_closure_ratio"}
    if not required.issubset(frame):
        raise ValueError("Observations require time_d and crack_closure_ratio")
    frame["split"] = _grouped_split(frame)
    evidence_label = (
        PUBLIC_CALIBRATED
        if "curation_status" in frame and len(frame)
        and frame["curation_status"].eq("approved").all()
        else UNCALIBRATED
    )
    empirical = _fit_empirical(frame.loc[frame["split"] == "train"])
    prediction_rows = []
    metric_rows = []
    observed = pd.to_numeric(frame["crack_closure_ratio"], errors="coerce").to_numpy(float)
    for structure in STRUCTURES:
        started = time.perf_counter()
        if structure == "first_order_empirical":
            times = pd.to_numeric(frame["time_d"], errors="coerce").to_numpy(float)
            predicted = (empirical[0] * (1.0 - np.exp(-empirical[1] * times))) if empirical is not None else np.full(len(frame), np.nan)
        elif structure == "zero_repair":
            predicted = np.zeros(len(frame))
        else:
            predicted = predictor(_structure_config(base, structure), frame)["crack_closure_ratio"].to_numpy(float)
        runtime = time.perf_counter() - started
        for index, value in enumerate(predicted):
            prediction_rows.append({
                "structure": structure, "row_index": index, "specimen_id": frame.get("specimen_id", pd.Series(index=frame.index, dtype=object)).iloc[index],
                "split": frame["split"].iloc[index], "time_d": frame["time_d"].iloc[index],
                "observed": observed[index], "predicted": value,
            })
        for split in ("train", "internal_test"):
            mask = frame["split"].to_numpy() == split
            metrics = information_criteria(observed[mask], predicted[mask], PARAMETER_COUNTS[structure])
            metric_rows.append({
                "structure": structure, "split": split, "parameter_count": PARAMETER_COUNTS[structure],
                "runtime_s": runtime, **metrics,
                "answers_design_question": structure not in {"first_order_empirical", "zero_repair"},
                "evidence_label": evidence_label,
            })
    pd.DataFrame(prediction_rows).to_csv(output_dir / "model_structure_predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(output_dir / "model_structure_comparison.csv", index=False)
    summary = {
        "structures": list(STRUCTURES), "observations_supplied": True,
        "information_criteria_available": True,
        "empirical_fit": empirical.tolist() if empirical is not None else None,
        "condition_resolved": True, "grouped_by_specimen": True,
    }
    (output_dir / "model_structure_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

"""Practical identifiability diagnostics distinct from global sensitivity."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .analysis import _get_parameter, _set_parameter
from .config import ModelConfig
from .model import simulate_0d


DEFAULT_PARAMETERS = (
    "kinetics.capsule_release_s",
    "kinetics.activity_multiplier",
    "kinetics.k_oxygen_mol_m3",
    "chemistry.calcite_rate_mol_m3_s",
    "chemistry.wall_deposition_fraction",
    "kinetics.capsule_csh_volume_fraction",
)
DEFAULT_OUTPUTS = (
    "lactate_mol_m3", "calcite_mol_m3", "crack_closure_ratio",
    "permeability_ratio", "csh_volume_fraction",
)
MEASUREMENT_MAP = {
    "kinetics.capsule_release_s": "released substrate concentration over the first 24 h",
    "kinetics.activity_multiplier": "substrate depletion or aggregate activity time series",
    "kinetics.k_oxygen_mol_m3": "paired dissolved oxygen and substrate depletion",
    "chemistry.calcite_rate_mol_m3_s": "time-resolved CaCO3 mass with calcium and DIC",
    "chemistry.wall_deposition_fraction": "paired CaCO3 mass and local crack-width microscopy",
    "kinetics.capsule_csh_volume_fraction": "C-S-H payload/release and abiotic closure control",
}


def matrix_diagnostics(
    sensitivity: np.ndarray, parameter_names: Sequence[str], rank_tolerance: float = 1e-8
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    """Diagnose a standardized local sensitivity matrix with an FIM/SVD."""

    matrix = np.asarray(sensitivity, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != len(parameter_names):
        raise ValueError("Sensitivity matrix columns must match parameter names")
    fim = matrix.T @ matrix
    eigenvalues = np.linalg.eigvalsh(fim)[::-1]
    singular = np.linalg.svd(matrix, compute_uv=False)
    threshold = (singular[0] * rank_tolerance) if len(singular) else 0.0
    rank = int(np.sum(singular > threshold))
    condition = float(singular[0] / singular[-1]) if len(singular) and singular[-1] > threshold else np.inf
    covariance = np.linalg.pinv(fim, rcond=rank_tolerance)
    scale = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(scale, scale)
    correlation = np.divide(covariance, denominator, out=np.zeros_like(covariance), where=denominator > 0)
    np.fill_diagonal(correlation, 1.0)
    norms = np.linalg.norm(matrix, axis=0)
    max_norm = max(float(norms.max()) if len(norms) else 0.0, 1e-30)
    rows = []
    for index, name in enumerate(parameter_names):
        relative = float(norms[index] / max_norm)
        max_correlation = float(np.max(np.abs(np.delete(correlation[index], index)))) if len(parameter_names) > 1 else 0.0
        if relative < 1e-4 or rank == 0:
            label = "not_estimable"
        elif relative < 0.05 or max_correlation >= 0.98 or rank < len(parameter_names):
            label = "weakly_estimable"
        else:
            label = "estimable"
        rows.append({
            "parameter": name, "sensitivity_norm": float(norms[index]),
            "relative_sensitivity": relative, "max_abs_correlation": max_correlation,
            "identifiability_label": label,
        })
    correlation_frame = pd.DataFrame(correlation, index=parameter_names, columns=parameter_names)
    summary = {
        "parameter_count": len(parameter_names), "fim_rank": rank,
        "full_rank": rank == len(parameter_names), "condition_number": condition,
        "fim_eigenvalues": eigenvalues.tolist(),
        "unidentifiable_combinations": max(len(parameter_names) - rank, 0),
    }
    return pd.DataFrame(rows), correlation_frame, summary


def _observable_vector(
    config: ModelConfig, times_d: Sequence[float], outputs: Sequence[str]
) -> Tuple[np.ndarray, Sequence[str]]:
    trial = copy.deepcopy(config)
    trial.simulation.days = max(float(max(times_d)), 0.01)
    trial.simulation.output_interval_days = min(max(trial.simulation.days / 80.0, 0.01), 1.0)
    frame = simulate_0d(trial).frame
    values, labels = [], []
    for output in outputs:
        for time_d in times_d:
            values.append(float(np.interp(time_d, frame["time_d"], frame[output])))
            labels.append("{}@{}d".format(output, time_d))
    return np.asarray(values), labels


def identifiability_analysis(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    parameters: Sequence[str] = DEFAULT_PARAMETERS,
    outputs: Sequence[str] = DEFAULT_OUTPUTS,
    times_d: Sequence[float] = (1.0, 7.0, 14.0, 28.0),
    relative_step: float = 0.02,
    measurement_scales: Optional[Mapping[str, float]] = None,
) -> Dict[str, object]:
    """Generate prior-design local/FIM diagnostics without claiming calibration."""

    base = copy.deepcopy(config or ModelConfig())
    baseline, labels = _observable_vector(base, times_d, outputs)
    scales = []
    for output in outputs:
        selected = baseline[len(scales):len(scales) + len(times_d)]
        default_scale = max(float(np.ptp(selected)), float(np.max(np.abs(selected))) * 0.05, 1e-8)
        scales.extend([float((measurement_scales or {}).get(output, default_scale))] * len(times_d))
    scales_array = np.asarray(scales)
    matrix = np.empty((len(baseline), len(parameters)))
    records = []
    for column, parameter in enumerate(parameters):
        value = _get_parameter(base, parameter)
        low, high = value * (1.0 - relative_step), value * (1.0 + relative_step)
        if value == 0:
            high, low = relative_step, -relative_step
        lower_config, upper_config = copy.deepcopy(base), copy.deepcopy(base)
        _set_parameter(lower_config, parameter, low)
        _set_parameter(upper_config, parameter, high)
        lower_values, _ = _observable_vector(lower_config, times_d, outputs)
        upper_values, _ = _observable_vector(upper_config, times_d, outputs)
        derivative_log = (upper_values - lower_values) / (2.0 * relative_step)
        matrix[:, column] = derivative_log / scales_array
        for row, label in enumerate(labels):
            records.append({
                "observation": label, "parameter": parameter,
                "standardized_log_sensitivity": matrix[row, column],
                "measurement_scale": scales_array[row],
                "evidence_class": "literature_prior_design_diagnostic",
            })
    table, correlation, summary = matrix_diagnostics(matrix, parameters)
    table["required_measurement"] = table["parameter"].map(MEASUREMENT_MAP).fillna("parameter-specific time series")
    measurement = table[["parameter", "identifiability_label", "required_measurement"]].copy()
    measurement["priority"] = measurement["identifiability_label"].map({
        "not_estimable": "highest", "weakly_estimable": "high", "estimable": "confirmatory",
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    pd.DataFrame(records).to_csv(output_dir / "local_sensitivity.csv", index=False)
    correlation.to_csv(output_dir / "parameter_correlation.csv")
    table.to_csv(output_dir / "identifiability_table.csv", index=False)
    measurement.to_csv(output_dir / "recommended_measurements.csv", index=False)
    summary.update({
        "evidence_class": "literature_prior_design_diagnostic",
        "sensitivity_is_not_identifiability": True,
        "profile_likelihood": "not_executable_without_observed_likelihood",
        "calibrated": False,
    })
    (output_dir / "identifiability_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary

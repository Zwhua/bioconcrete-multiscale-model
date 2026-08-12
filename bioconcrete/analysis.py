"""Calibration and global sensitivity analysis for model parameters."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import qmc

from .config import ModelConfig, PARAMETER_PROVENANCE
from .model import simulate_0d


CALIBRATION_COLUMNS = (
    "time_d",
    "group",
    "crack_width_mm",
    "pH",
    "oxygen_mg_L",
    "lactate_mM",
    "cfu_mL",
    "caco3_mg",
    "healing_ratio",
    "permeability_ratio",
    "sorptivity",
)


def _get_parameter(config: ModelConfig, path: str) -> float:
    section, name = path.split(".", 1)
    return float(getattr(getattr(config, section), name))


def _set_parameter(config: ModelConfig, path: str, value: float) -> None:
    section, name = path.split(".", 1)
    setattr(getattr(config, section), name, float(value))


def validate_experiment_csv(frame: pd.DataFrame) -> None:
    if "time_d" not in frame.columns:
        raise ValueError("Experiment CSV must contain time_d")
    usable = {"lactate_mM", "cfu_mL", "caco3_mg", "healing_ratio", "permeability_ratio", "sorptivity"}
    if not (usable & set(frame.columns)):
        raise ValueError("Experiment CSV contains no supported response column")
    if (pd.to_numeric(frame["time_d"], errors="coerce") < 0).any():
        raise ValueError("time_d cannot be negative")


def _model_predictions(config: ModelConfig, observations: pd.DataFrame) -> Dict[str, np.ndarray]:
    config.simulation.days = max(float(observations["time_d"].max()), 0.01)
    config.simulation.output_interval_days = max(config.simulation.days / 100.0, 0.01)
    model = simulate_0d(config).frame
    times = observations["time_d"].to_numpy(dtype=float)
    interpolate = lambda column: np.interp(times, model["time_d"], model[column])
    crack_volume_m3 = config.transport.crack_width_mm / 1000.0
    return {
        "lactate_mM": interpolate("lactate_mol_m3"),
        "cfu_mL": interpolate("active_density_rel"),
        "caco3_mg": interpolate("calcite_mol_m3") * config.chemistry.calcite_molar_mass_kg_mol * crack_volume_m3 * 1e6,
        "healing_ratio": interpolate("healing_ratio"),
        "permeability_ratio": interpolate("permeability_ratio"),
        "sorptivity": interpolate("sorptivity_ratio"),
    }


def calibrate(
    data_path: Path,
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    parameters: Sequence[str] = (
        "kinetics.qmax_lactate_mol_m3_s",
        "kinetics.capsule_release_s",
        "chemistry.calcite_rate_mol_m3_s",
    ),
    bootstrap_samples: int = 20,
) -> Dict[str, object]:
    frame = pd.read_csv(data_path)
    validate_experiment_csv(frame)
    base = copy.deepcopy(config or ModelConfig())
    bounds = [PARAMETER_PROVENANCE[name][2:4] for name in parameters]
    lower = np.log([item[0] for item in bounds])
    upper = np.log([item[1] for item in bounds])
    initial = np.log([_get_parameter(base, name) for name in parameters])
    response_columns = [
        name for name in ("lactate_mM", "cfu_mL", "caco3_mg", "healing_ratio", "permeability_ratio", "sorptivity")
        if name in frame.columns and frame[name].notna().any()
    ]

    scales = {}
    for column in response_columns:
        values = pd.to_numeric(frame[column], errors="coerce")
        scales[column] = max(float(values.std()), float(values.abs().median()) * 0.1, 1e-6)

    def residual(log_values: np.ndarray, sample: pd.DataFrame) -> np.ndarray:
        trial = copy.deepcopy(base)
        for name, value in zip(parameters, np.exp(log_values)):
            _set_parameter(trial, name, value)
        predicted = _model_predictions(trial, sample)
        residuals = []
        for column in response_columns:
            observed = pd.to_numeric(sample[column], errors="coerce").to_numpy(dtype=float)
            mask = np.isfinite(observed)
            residuals.extend(((predicted[column][mask] - observed[mask]) / scales[column]).tolist())
        return np.asarray(residuals)

    fit = least_squares(residual, initial, bounds=(lower, upper), args=(frame,), max_nfev=80)
    fitted_values = np.exp(fit.x)
    rng = np.random.RandomState(base.simulation.random_seed)
    bootstrap = []
    for _ in range(max(bootstrap_samples, 0)):
        sampled = frame.iloc[rng.randint(0, len(frame), len(frame))].sort_values("time_d").reset_index(drop=True)
        try:
            result = least_squares(residual, fit.x, bounds=(lower, upper), args=(sampled,), max_nfev=30)
            bootstrap.append(np.exp(result.x))
        except (RuntimeError, ValueError):
            continue
    intervals = np.asarray(bootstrap) if bootstrap else np.empty((0, len(parameters)))
    rows = []
    for index, (name, value) in enumerate(zip(parameters, fitted_values)):
        rows.append(
            {
                "parameter": name,
                "estimate": value,
                "ci_low": float(np.quantile(intervals[:, index], 0.025)) if len(intervals) else np.nan,
                "ci_high": float(np.quantile(intervals[:, index], 0.975)) if len(intervals) else np.nan,
            }
        )
        _set_parameter(base, name, value)
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "fitted_parameters.csv", index=False)
    base.save(output_dir / "fitted_config.json")
    summary = {
        "success": bool(fit.success),
        "message": fit.message,
        "cost": float(fit.cost),
        "observations": int(len(frame)),
        "responses": response_columns,
        "bootstrap_successes": int(len(intervals)),
    }
    (output_dir / "calibration_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def _sample_parameter(unit_value: float, lower: float, upper: float) -> float:
    if lower > 0 and upper / lower > 20:
        return float(np.exp(np.log(lower) + unit_value * (np.log(upper) - np.log(lower))))
    return float(lower + unit_value * (upper - lower))


def _evaluate(config: ModelConfig, paths: Sequence[str], unit_values: np.ndarray) -> float:
    trial = copy.deepcopy(config)
    for path, unit_value in zip(paths, unit_values):
        lower, upper = PARAMETER_PROVENANCE[path][2:4]
        _set_parameter(trial, path, _sample_parameter(float(unit_value), lower, upper))
    return simulate_0d(trial).summary["mean_healing_ratio"]


def sensitivity(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    samples: int = 8,
    parameters: Sequence[str] = (
        "kinetics.qmax_lactate_mol_m3_s",
        "kinetics.k_lactate_mol_m3",
        "kinetics.k_oxygen_mol_m3",
        "kinetics.decay_s",
        "kinetics.capsule_release_s",
        "chemistry.calcite_rate_mol_m3_s",
        "chemistry.portlandite_mol_m3",
        "transport.crack_width_mm",
    ),
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    config = copy.deepcopy(config or ModelConfig())
    dimension = len(parameters)
    samples = max(int(samples), 4)
    sampler = qmc.Sobol(d=2 * dimension, scramble=True, seed=config.simulation.random_seed)
    power = int(np.ceil(np.log2(samples)))
    points = sampler.random_base2(power)[:samples]
    matrix_a, matrix_b = points[:, :dimension], points[:, dimension:]
    y_a = np.array([_evaluate(config, parameters, row) for row in matrix_a])
    y_b = np.array([_evaluate(config, parameters, row) for row in matrix_b])
    variance = max(float(np.var(np.concatenate([y_a, y_b]), ddof=1)), 1e-30)
    sobol_rows = []
    morris_rows = []
    delta = 0.1
    for index, parameter in enumerate(parameters):
        mixed_outputs = []
        effects = []
        for row_index in range(samples):
            mixed = matrix_a[row_index].copy()
            mixed[index] = matrix_b[row_index, index]
            y_mixed = _evaluate(config, parameters, mixed)
            mixed_outputs.append(y_mixed)
            denominator = matrix_b[row_index, index] - matrix_a[row_index, index]
            if abs(denominator) > 1e-8:
                effects.append((y_mixed - y_a[row_index]) / denominator)
        y_mixed_array = np.asarray(mixed_outputs)
        first_order = float(np.mean(y_b * (y_mixed_array - y_a)) / variance)
        total_order = float(0.5 * np.mean((y_a - y_mixed_array) ** 2) / variance)
        sobol_rows.append({"parameter": parameter, "S1": first_order, "ST": total_order})
        morris_rows.append(
            {
                "parameter": parameter,
                "mu": float(np.mean(effects)) if effects else np.nan,
                "mu_star": float(np.mean(np.abs(effects))) if effects else np.nan,
                "sigma": float(np.std(effects, ddof=1)) if len(effects) > 1 else np.nan,
            }
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    sobol = pd.DataFrame(sobol_rows)
    morris = pd.DataFrame(morris_rows)
    sobol.to_csv(output_dir / "sobol_indices.csv", index=False)
    morris.to_csv(output_dir / "morris_indices.csv", index=False)
    return morris, sobol

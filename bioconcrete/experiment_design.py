"""Numerical model-driven greedy D-optimal experimental design."""

from __future__ import annotations

import copy
import hashlib
import itertools
import json
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .analysis import _get_parameter, _set_parameter
from .config import ModelConfig
from .evidence_state import UNCALIBRATED
from .model import simulate_0d


DESIGN_PARAMETERS = (
    "kinetics.capsule_release_s",
    "kinetics.activity_multiplier",
    "chemistry.calcite_rate_mol_m3_s",
    "chemistry.wall_deposition_fraction",
    "kinetics.capsule_csh_volume_fraction",
)
OBSERVABLES = {
    "crack_closure_ratio": ("crack_closure_ratio", 0.02),
    "calcite_mass_mg": ("calcite_mol_m3", 0.01),
    "substrate_mol_m3": ("lactate_mol_m3", 50.0),
    "oxygen_mol_m3": ("oxygen_mol_m3", 0.01),
    "calcium_mol_m3": ("calcium_mol_m3", 10.0),
    "ph": ("ph", 0.10),
    "csh_fill_fraction": ("csh_volume_fraction", 0.002),
}
SYSTEMS = ("complete", "no_csh", "abiotic_csh", "no_biological_mineralization")


def d_optimal_score(sensitivity: np.ndarray, prior_information: Optional[np.ndarray] = None) -> float:
    """Return log determinant for a sensitivity matrix and prior information."""

    matrix = np.atleast_2d(np.asarray(sensitivity, dtype=float))
    prior = np.eye(matrix.shape[1]) * 1e-6 if prior_information is None else np.asarray(prior_information)
    sign, value = np.linalg.slogdet(prior + matrix.T @ matrix)
    return float(value) if sign > 0 else -np.inf


def fisher_information(jacobian: np.ndarray, covariance: np.ndarray) -> np.ndarray:
    """Compute J.T Sigma^-1 J with explicit measurement covariance."""

    jacobian = np.atleast_2d(np.asarray(jacobian, dtype=float))
    covariance = np.atleast_2d(np.asarray(covariance, dtype=float))
    return jacobian.T @ np.linalg.pinv(covariance) @ jacobian


def _apply_system(config: ModelConfig, system: str) -> None:
    if system == "no_csh":
        config.kinetics.capsule_csh_volume_fraction = 0.0
    elif system == "abiotic_csh":
        config.kinetics.spore_density_rel = 0.0
        config.kinetics.active_density_rel = 0.0
        config.chemistry.calcite_rate_mol_m3_s = 0.0
    elif system == "no_biological_mineralization":
        config.kinetics.spore_density_rel = 0.0
        config.kinetics.active_density_rel = 0.0
        config.chemistry.calcite_rate_mol_m3_s = 0.0
        config.kinetics.capsule_csh_volume_fraction = 0.0
    elif system != "complete":
        raise ValueError("Unknown experimental system: {}".format(system))


def _configure(base: ModelConfig, width: float, dose: float, wet: float, system: str) -> ModelConfig:
    trial = copy.deepcopy(base)
    trial.transport.crack_width_mm = width
    trial.kinetics.agent_dosage_multiplier = dose
    trial.environment.wet_hours_per_day = wet
    trial.simulation.days = 28.0
    trial.simulation.output_interval_days = 1.0
    _apply_system(trial, system)
    trial.validate()
    return trial


def _trajectory(config: ModelConfig) -> pd.DataFrame:
    return simulate_0d(config).frame


def _value(frame: pd.DataFrame, observable: str, time_d: float, config: ModelConfig) -> float:
    column, _ = OBSERVABLES[observable]
    value = float(np.interp(time_d, frame["time_d"], frame[column]))
    if observable == "calcite_mass_mg":
        volume = (
            config.transport.crack_length_mm * config.transport.crack_width_mm
            * config.transport.crack_depth_mm * 1e-9
        )
        value *= config.chemistry.calcite_molar_mass_kg_mol * volume * 1e6
    return value


def numerical_candidate_table(
    config: Optional[ModelConfig] = None,
    relative_step: float = 0.05,
    trajectory_runner: Callable[[ModelConfig], pd.DataFrame] = _trajectory,
    widths: Sequence[float] = (0.1, 0.3, 0.5),
    doses: Sequence[float] = (0.5, 1.0, 2.0),
    wettings: Sequence[float] = (6.0, 12.0, 24.0),
    times: Sequence[float] = (1.0, 7.0, 14.0, 21.0, 28.0),
    systems: Sequence[str] = SYSTEMS,
) -> pd.DataFrame:
    """Generate candidates from cached model trajectories and central differences."""

    if relative_step <= 0 or relative_step >= 1:
        raise ValueError("relative_step must be in (0, 1)")
    base = copy.deepcopy(config or ModelConfig())
    rows: List[Dict[str, object]] = []
    for width, dose, wet, system in itertools.product(widths, doses, wettings, systems):
        condition = _configure(base, width, dose, wet, system)
        baseline = trajectory_runner(condition)
        derivative_frames: Dict[str, Tuple[pd.DataFrame, pd.DataFrame, ModelConfig, ModelConfig]] = {}
        for parameter in DESIGN_PARAMETERS:
            value = _get_parameter(condition, parameter)
            minus, plus = copy.deepcopy(condition), copy.deepcopy(condition)
            _set_parameter(minus, parameter, value * (1.0 - relative_step))
            _set_parameter(plus, parameter, value * (1.0 + relative_step))
            minus.validate()
            plus.validate()
            derivative_frames[parameter] = (
                trajectory_runner(minus), trajectory_runner(plus), minus, plus
            )
        for time_d, observable in itertools.product(times, OBSERVABLES):
            _, sigma = OBSERVABLES[observable]
            vector = []
            for parameter in DESIGN_PARAMETERS:
                minus_frame, plus_frame, minus_config, plus_config = derivative_frames[parameter]
                low = _value(minus_frame, observable, time_d, minus_config)
                high = _value(plus_frame, observable, time_d, plus_config)
                vector.append((high - low) / (2.0 * relative_step))
            identifier = hashlib.sha256(
                "{}|{}|{}|{}|{}|{}".format(width, dose, wet, system, time_d, observable).encode()
            ).hexdigest()[:16]
            execution_cost = 1.0 + time_d / 28.0 + (0.5 if system != "complete" else 0.0)
            rows.append({
                "experiment_id": identifier, "crack_width_mm": width,
                "agent_dosage": dose, "wet_hours_per_day": wet,
                "system": system, "time_d": time_d, "observable": observable,
                "predicted_value": _value(baseline, observable, time_d, condition),
                "measurement_sd": sigma, "measurement_error_source": "literature_prior_pending_public_measurement_fit",
                "execution_cost": execution_cost, "biological_replicates": 3,
                **{"sensitivity_{}".format(name): value for name, value in zip(DESIGN_PARAMETERS, vector)},
                "evidence_label": UNCALIBRATED,
            })
    return pd.DataFrame(rows)


def greedy_d_optimal(candidates: pd.DataFrame, count: int = 10) -> Tuple[pd.DataFrame, np.ndarray]:
    """Select a complementary experiment set by sequential Fisher information gain."""

    sensitivity_columns = ["sensitivity_{}".format(name) for name in DESIGN_PARAMETERS]
    information = np.eye(len(DESIGN_PARAMETERS)) * 1e-6
    selected: List[int] = []
    available = list(candidates.index)
    records = []
    for rank in range(1, min(count, len(available)) + 1):
        current_score = d_optimal_score(np.empty((0, len(DESIGN_PARAMETERS))), information)
        best = None
        for index in available:
            row = candidates.loc[index]
            jacobian = row[sensitivity_columns].to_numpy(float).reshape(1, -1)
            covariance = np.array([[float(row["measurement_sd"]) ** 2]])
            candidate_information = information + fisher_information(jacobian, covariance)
            gain = d_optimal_score(np.empty((0, len(DESIGN_PARAMETERS))), candidate_information) - current_score
            utility = gain / float(row["execution_cost"])
            if best is None or utility > best[0]:
                best = (utility, gain, index, candidate_information)
        assert best is not None
        utility, gain, index, information = best
        available.remove(index)
        selected.append(index)
        records.append({"rank": rank, "information_gain": gain, "cost_adjusted_gain": utility})
    result = candidates.loc[selected].reset_index(drop=True).join(pd.DataFrame(records))
    return result, information


def rank_experiments(
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    trajectory_runner: Callable[[ModelConfig], pd.DataFrame] = _trajectory,
    smoke: bool = False,
) -> Dict[str, object]:
    """Run numerical sensitivities and select minimum and ideal D-optimal plans."""

    options = ({"widths": (0.3,), "doses": (1.0,), "wettings": (12.0,),
                "times": (1.0, 28.0), "systems": ("complete",)} if smoke else {})
    frame = numerical_candidate_table(config, trajectory_runner=trajectory_runner, **options)
    selected, final_information = greedy_d_optimal(frame, 10)
    selected["plan"] = np.where(selected["rank"] <= 5, "minimum_executable", "ideal_extension")
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "candidate_experiments.csv", index=False)
    selected.to_csv(output_dir / "recommended_experiments.csv", index=False)
    np.savetxt(output_dir / "selected_fisher_information.csv", final_information, delimiter=",")
    summary = {
        "method": "numerical greedy D-optimal",
        "candidate_count": len(frame), "minimum_executable_count": 5,
        "ideal_plan_count": 10, "parameters": list(DESIGN_PARAMETERS),
        "measurement_covariance": "diagonal literature prior; public measurement fit pending",
        "evidence_label": UNCALIBRATED, "experimental_validation": False,
        "smoke_test": smoke,
    }
    (output_dir / "experiment_design_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

"""Calibration, frozen external validation, and measurement uncertainty."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .analysis import _get_parameter, _set_parameter
from .config import ModelConfig, PARAMETER_PROVENANCE
from .model import simulate_0d


SHARED_PARAMETERS = (
    "kinetics.capsule_release_s",
    "kinetics.activity_multiplier",
    "chemistry.calcite_rate_mol_m3_s",
    "chemistry.wall_deposition_fraction",
)


def _digest_config(config: ModelConfig) -> str:
    payload = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _metrics(observed: np.ndarray, predicted: np.ndarray) -> Dict[str, float]:
    mask = np.isfinite(observed) & np.isfinite(predicted)
    if not mask.any():
        return {"n": 0, "mae": np.nan, "rmse": np.nan, "r2": np.nan}
    obs, pred = observed[mask], predicted[mask]
    residual = pred - obs
    denominator = float(np.sum((obs - np.mean(obs)) ** 2))
    result = {
        "n": int(mask.sum()),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0 else np.nan,
    }
    rss = float(np.sum(residual**2))
    result["aic"] = float(mask.sum() * np.log(max(rss / mask.sum(), 1e-30)))
    return result


def _closure(frame: pd.DataFrame) -> np.ndarray:
    initial = pd.to_numeric(frame["initial_crack_width_mm"], errors="coerce").to_numpy(float)
    current = pd.to_numeric(frame["current_crack_width_mm"], errors="coerce").to_numpy(float)
    return np.clip(1.0 - current / initial, 0.0, 1.0)


def _predict_rows(config: ModelConfig, observations: pd.DataFrame) -> np.ndarray:
    observations = observations.reset_index(drop=True).copy()
    for column in ("initial_crack_width_mm", "wet_hours_per_day", "agent_dosage", "time_d"):
        observations[column] = pd.to_numeric(observations[column], errors="coerce")
    predictions = np.full(len(observations), np.nan)
    for (width, wet, dosage), indexes in observations.groupby(
        ["initial_crack_width_mm", "wet_hours_per_day", "agent_dosage"], dropna=False
    ).groups.items():
        trial = copy.deepcopy(config)
        if np.isfinite(width):
            trial.transport.crack_width_mm = float(width)
        if np.isfinite(wet):
            trial.environment.wet_hours_per_day = float(wet)
        if np.isfinite(dosage):
            trial.kinetics.capsule_calcium_lactate_mol_m3 *= max(float(dosage), 0.0)
        subset = observations.loc[indexes]
        days = max(float(pd.to_numeric(subset["time_d"], errors="coerce").max()), 0.01)
        trial.simulation.days = days
        trial.simulation.output_interval_days = max(days / 80.0, 0.05)
        model = simulate_0d(trial).frame
        predictions[np.asarray(list(indexes), dtype=int)] = np.interp(
            pd.to_numeric(subset["time_d"], errors="coerce"),
            model["time_d"], model["crack_closure_ratio"],
        )
    return predictions


def calibrate_public(data_path: Path, output_dir: Path, config: Optional[ModelConfig] = None,
                     bootstrap_samples: int = 20, profile_points: int = 0) -> Dict[str, object]:
    frame = pd.read_csv(data_path)
    required = {"dataset_id", "specimen_id", "split", "time_d", "initial_crack_width_mm", "current_crack_width_mm"}
    if not required.issubset(frame.columns):
        raise ValueError("Public observation table is missing required columns")
    train = frame.loc[frame["split"] == "train"].copy()
    holdout = frame.loc[frame["split"] == "internal_test"].copy()
    if set(train["specimen_id"]) & set(holdout["specimen_id"]):
        raise ValueError("Specimen leakage detected between train and holdout")
    observed = _closure(train)
    usable = np.isfinite(observed) & np.isfinite(pd.to_numeric(train["time_d"], errors="coerce"))
    train = train.loc[usable].reset_index(drop=True)
    observed = observed[usable]
    base = copy.deepcopy(config or ModelConfig())
    lower, upper, initial = [], [], []
    for name in SHARED_PARAMETERS:
        if name in PARAMETER_PROVENANCE:
            bounds = PARAMETER_PROVENANCE[name][2:4]
        else:
            value = _get_parameter(base, name)
            bounds = (max(value / 10.0, 1e-8), value * 10.0)
        lower.append(bounds[0]); upper.append(bounds[1]); initial.append(_get_parameter(base, name))

    def residual(log_values: np.ndarray, sample: pd.DataFrame, target: np.ndarray) -> np.ndarray:
        trial = copy.deepcopy(base)
        for name, value in zip(SHARED_PARAMETERS, np.exp(log_values)):
            _set_parameter(trial, name, value)
        return _predict_rows(trial, sample) - target

    active_indexes = [0, 1, 3]
    calcite_values = pd.to_numeric(train.get("calcite_mass_mg"), errors="coerce")
    calcite_observed = bool(calcite_values.notna().sum() >= 2)
    if calcite_observed:
        active_indexes.insert(2, 2)
    if len(train) < 2 * len(active_indexes):
        raise ValueError("Each fitted parameter requires at least two independent observations")
    current = np.log(initial)
    stages = []
    early_limit = float(pd.to_numeric(train["time_d"], errors="coerce").median())
    stage_specs = [
        ("release_and_effective_activity", [0, 1], train.loc[train["time_d"] <= early_limit].reset_index(drop=True)),
        ("precipitation", [2] if calcite_observed else [], train),
        ("wall_geometry_response", [3], train),
    ]
    for stage_name, indexes, sample in stage_specs:
        indexes = [index for index in indexes if index in active_indexes]
        if not indexes or len(sample) < 2 * len(indexes):
            stages.append({"stage": stage_name, "status": "fixed_to_prior", "parameters": [SHARED_PARAMETERS[i] for i in indexes]})
            continue
        target = _closure(sample)
        def stage_residual(stage_values: np.ndarray) -> np.ndarray:
            values = current.copy()
            values[indexes] = stage_values
            return residual(values, sample, target)
        stage_fit = least_squares(
            stage_residual, current[indexes],
            bounds=(np.asarray(np.log(lower))[indexes], np.asarray(np.log(upper))[indexes]), max_nfev=80,
        )
        current[indexes] = stage_fit.x
        stages.append({"stage": stage_name, "status": "fitted", "parameters": [SHARED_PARAMETERS[i] for i in indexes],
                       "cost": float(stage_fit.cost), "success": bool(stage_fit.success)})
    fit = least_squares(
        lambda values: residual(
            np.asarray([values[active_indexes.index(i)] if i in active_indexes else current[i]
                        for i in range(len(SHARED_PARAMETERS))]), train, observed
        ),
        current[active_indexes],
        bounds=(np.asarray(np.log(lower))[active_indexes], np.asarray(np.log(upper))[active_indexes]),
        max_nfev=60,
    )
    final_values = current.copy()
    final_values[active_indexes] = fit.x
    fit.x = final_values
    fitted = copy.deepcopy(base)
    for name, value in zip(SHARED_PARAMETERS, np.exp(fit.x)):
        _set_parameter(fitted, name, value)
    rng = np.random.RandomState(base.simulation.random_seed)
    boot = []
    specimens = train["specimen_id"].astype(str).unique()
    for _ in range(max(bootstrap_samples, 0)):
        sampled_ids = rng.choice(specimens, len(specimens), replace=True)
        sample = pd.concat([train.loc[train["specimen_id"].astype(str) == item] for item in sampled_ids], ignore_index=True)
        target = _closure(sample)
        try:
            def bootstrap_residual(active_values: np.ndarray) -> np.ndarray:
                values = fit.x.copy()
                values[active_indexes] = active_values
                return residual(values, sample, target)
            result = least_squares(
                bootstrap_residual, fit.x[active_indexes],
                bounds=(np.asarray(np.log(lower))[active_indexes], np.asarray(np.log(upper))[active_indexes]),
                max_nfev=30,
            )
            values = fit.x.copy()
            values[active_indexes] = result.x
            boot.append(np.exp(values))
        except (RuntimeError, ValueError):
            pass
    output_dir.mkdir(parents=True, exist_ok=True)
    fitted.save(output_dir / "frozen_config.json")
    train_predictions = _predict_rows(fitted, train)
    holdout_observed = _closure(holdout) if len(holdout) else np.array([])
    holdout_predictions = _predict_rows(fitted, holdout) if len(holdout) else np.array([])
    intervals = np.asarray(boot)
    parameters = []
    for index, name in enumerate(SHARED_PARAMETERS):
        parameters.append({
            "parameter": name, "estimate": _get_parameter(fitted, name),
            "ci_low": float(np.quantile(intervals[:, index], .025)) if len(intervals) else np.nan,
            "ci_high": float(np.quantile(intervals[:, index], .975)) if len(intervals) else np.nan,
        })
    pd.DataFrame(parameters).to_csv(output_dir / "fitted_parameters.csv", index=False)
    if len(intervals):
        bootstrap_frame = pd.DataFrame(intervals, columns=SHARED_PARAMETERS)
        bootstrap_frame.to_csv(output_dir / "bootstrap_parameters.csv", index=False)
        bootstrap_frame.corr().to_csv(output_dir / "parameter_correlation.csv")
        prediction_rows = []
        prediction_days = np.array([7.0, 14.0, 21.0, 28.0])
        ensemble_predictions = []
        for values in intervals[:200]:
            trial = copy.deepcopy(fitted)
            for name, value in zip(SHARED_PARAMETERS, values):
                _set_parameter(trial, name, float(value))
            trial.simulation.days = 28.0
            trial.simulation.output_interval_days = 1.0
            result = simulate_0d(trial).frame
            ensemble_predictions.append(np.interp(
                prediction_days, result["time_d"], result["crack_closure_ratio"]
            ))
        ensemble_predictions = np.asarray(ensemble_predictions)
        for index, day in enumerate(prediction_days):
            prediction_rows.append({
                "time_d": day,
                "closure_median": float(np.quantile(ensemble_predictions[:, index], .5)),
                "closure_low_95": float(np.quantile(ensemble_predictions[:, index], .025)),
                "closure_high_95": float(np.quantile(ensemble_predictions[:, index], .975)),
                "evidence_class": "model_prediction",
            })
        pd.DataFrame(prediction_rows).to_csv(output_dir / "calibrated_prediction_intervals.csv", index=False)
    profile_rows = []
    if profile_points >= 3:
        for fixed_index, fixed_name in enumerate(SHARED_PARAMETERS):
            grid = np.linspace(np.log(lower[fixed_index]), np.log(upper[fixed_index]), profile_points)
            free_indexes = [index for index in range(len(SHARED_PARAMETERS)) if index != fixed_index]
            for fixed_value in grid:
                def profiled_residual(free_values: np.ndarray) -> np.ndarray:
                    values = fit.x.copy()
                    values[fixed_index] = fixed_value
                    values[free_indexes] = free_values
                    return residual(values, train, observed)
                profile_fit = least_squares(
                    profiled_residual, fit.x[free_indexes],
                    bounds=(np.asarray(np.log(lower))[free_indexes], np.asarray(np.log(upper))[free_indexes]),
                    max_nfev=40,
                )
                profile_rows.append({"parameter": fixed_name, "fixed_value": float(np.exp(fixed_value)),
                                     "cost": float(profile_fit.cost), "success": bool(profile_fit.success)})
        pd.DataFrame(profile_rows).to_csv(output_dir / "profile_likelihood.csv", index=False)
    artifact = {
        "evidence_class": "public_calibration_data",
        "project_experiment_rows": 0,
        "dataset_ids": sorted(frame["dataset_id"].dropna().unique().tolist()),
        "frozen_config_sha256": _digest_config(fitted),
        "shared_parameters": list(SHARED_PARAMETERS),
        "train_metrics": _metrics(observed, train_predictions),
        "internal_test_metrics": _metrics(holdout_observed, holdout_predictions),
        "bootstrap_successes": len(boot),
        "identifiability_rule": "at least two observations per fitted parameter; inspect bootstrap CI and correlation",
        "profile_likelihood_points": int(profile_points),
        "calibration_stages": stages,
        "fixed_to_prior": [SHARED_PARAMETERS[2]] if not calcite_observed else [],
    }
    (output_dir / "frozen_run.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    return artifact


def validate_external(data_path: Path, frozen_run_dir: Path, output_dir: Path) -> Dict[str, object]:
    frame = pd.read_csv(data_path)
    if not (frame["split"] == "external_validation").all():
        raise ValueError("External data must be labelled external_validation")
    config = ModelConfig.load(frozen_run_dir / "frozen_config.json")
    artifact = json.loads((frozen_run_dir / "frozen_run.json").read_text(encoding="utf-8"))
    if _digest_config(config) != artifact["frozen_config_sha256"]:
        raise ValueError("Frozen configuration checksum mismatch")
    observed = _closure(frame)
    mechanism = _predict_rows(config, frame)
    bootstrap_path = frozen_run_dir / "bootstrap_parameters.csv"
    ensemble = []
    if bootstrap_path.exists():
        samples = pd.read_csv(bootstrap_path)
        for _, values in samples.head(200).iterrows():
            trial = copy.deepcopy(config)
            for name in SHARED_PARAMETERS:
                _set_parameter(trial, name, float(values[name]))
            ensemble.append(_predict_rows(trial, frame))
    if ensemble:
        ensemble_values = np.asarray(ensemble)
        prediction_low = np.quantile(ensemble_values, .025, axis=0)
        prediction_high = np.quantile(ensemble_values, .975, axis=0)
        finite = np.isfinite(observed) & np.isfinite(prediction_low) & np.isfinite(prediction_high)
        coverage = float(np.mean((observed[finite] >= prediction_low[finite]) &
                                 (observed[finite] <= prediction_high[finite]))) if finite.any() else np.nan
    else:
        prediction_low = np.full(len(frame), np.nan)
        prediction_high = np.full(len(frame), np.nan)
        coverage = np.nan
    times = pd.to_numeric(frame["time_d"], errors="coerce").to_numpy(float)
    zero = np.zeros(len(frame))
    train_rmse = artifact["train_metrics"]["rmse"]
    rate = max(-np.log(max(1.0 - min(train_rmse if np.isfinite(train_rmse) else .1, .95), .05)) / 28.0, 1e-6)
    empirical = 1.0 - np.exp(-rate * times)
    results = {
        "evidence_class": "external_validation_data",
        "parameters_refitted": False,
        "frozen_config_sha256": artifact["frozen_config_sha256"],
        "mechanistic": _metrics(observed, mechanism),
        "prediction_interval_coverage": coverage,
        "prediction_interval_samples": len(ensemble),
        "zero_mineralization_baseline": _metrics(observed, zero),
        "first_order_baseline": _metrics(observed, empirical),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"observed_closure": observed, "mechanistic_prediction": mechanism,
                  "prediction_low_95": prediction_low, "prediction_high_95": prediction_high,
                  "zero_baseline": zero, "first_order_baseline": empirical}).to_csv(output_dir / "predictions.csv", index=False)
    (output_dir / "external_validation.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results


def fit_measurement_error(data_path: Path, output_dir: Path) -> Dict[str, object]:
    frame = pd.read_csv(data_path)
    numeric = frame.select_dtypes(include=[np.number])
    reference = next((name for name in numeric if "reference" in name.lower() or "manual" in name.lower()), None)
    prediction = next((name for name in numeric if name != reference and ("predict" in name.lower() or "model" in name.lower())), None)
    if reference is None or prediction is None:
        raise ValueError("Measurement table needs identifiable reference/manual and prediction columns")
    residual = numeric[prediction] - numeric[reference]
    result = {"evidence_class": "measurement_error_data", "rows": int(residual.notna().sum()),
              "bias": float(residual.mean()), "measurement_sd": float(residual.std(ddof=1)),
              "kinetic_parameters_fitted": False}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "measurement_error.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result

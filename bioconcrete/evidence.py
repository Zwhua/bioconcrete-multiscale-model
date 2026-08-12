"""Public-data calibration, frozen validation, and measurement uncertainty.

Scientific calibration in this module is deliberately separate from synthetic
parameter-recovery tests. Each response is compared in its own physical unit
and standardized by an auditable uncertainty estimate before residuals are
combined.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

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

OUTPUT_ALIASES = {
    "lactate_mM": ("lactate_mM",),
    "calcite_mass_mg": ("calcite_mass_mg", "caco3_mg"),
    "crack_closure_ratio": ("crack_closure_ratio", "healing_ratio"),
    "permeability_ratio": ("permeability_ratio",),
    "ph": ("pH", "ph"),
    "activation_state": ("activation_state", "effective_activity"),
    "cumulative_activity_h": ("cumulative_activity_h",),
}

PRESET_SCALES = {
    "lactate_mM": 10.0,
    "calcite_mass_mg": 0.10,
    "crack_closure_ratio": 0.05,
    "permeability_ratio": 0.05,
    "ph": 0.20,
    "activation_state": 0.10,
    "cumulative_activity_h": 1.0,
}

STAGES = (
    ("release_and_effective_activity", (0, 1),
     ("lactate_mM", "activation_state", "cumulative_activity_h")),
    ("precipitation", (2,), ("calcite_mass_mg",)),
    ("wall_geometry_response", (3,), ("crack_closure_ratio",)),
    ("material_performance", (), ("permeability_ratio",)),
)

OUTPUT_PARAMETER_INDEXES = {
    "lactate_mM": (0, 1),
    "activation_state": (0, 1),
    "cumulative_activity_h": (0, 1),
    "calcite_mass_mg": (2,),
    "crack_closure_ratio": (3,),
    "permeability_ratio": (),
    "ph": (),
}


def _digest_config(config: ModelConfig) -> str:
    payload = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _metrics(
    observed: np.ndarray,
    predicted: np.ndarray,
    parameter_count: int = 0,
    prediction_low: Optional[np.ndarray] = None,
    prediction_high: Optional[np.ndarray] = None,
) -> Dict[str, float]:
    """Return metrics in the response's physical units."""

    obs = np.asarray(observed, dtype=float)
    pred = np.asarray(predicted, dtype=float)
    mask = np.isfinite(obs) & np.isfinite(pred)
    if not mask.any():
        return {
            "n": 0, "mae": np.nan, "rmse": np.nan, "r2": np.nan,
            "calibration_slope": np.nan, "calibration_intercept": np.nan,
            "rss": np.nan, "aic": np.nan, "aicc": np.nan,
            "prediction_interval_coverage": np.nan,
            "mean_prediction_interval_width": np.nan,
        }
    obs, pred = obs[mask], pred[mask]
    residual = pred - obs
    rss = float(np.sum(residual**2))
    n = int(mask.sum())
    denominator = float(np.sum((obs - np.mean(obs)) ** 2))
    k = max(int(parameter_count), 0)
    aic = float(n * np.log(max(rss / n, 1e-30)) + 2 * k)
    aicc = float(aic + 2 * k * (k + 1) / (n - k - 1)) if n > k + 1 else np.nan
    if n >= 2 and float(np.ptp(pred)) > 1e-15:
        slope, intercept = np.polyfit(pred, obs, 1)
    else:
        slope, intercept = np.nan, np.nan
    coverage, interval_width = np.nan, np.nan
    if prediction_low is not None and prediction_high is not None:
        low = np.asarray(prediction_low, dtype=float)[mask]
        high = np.asarray(prediction_high, dtype=float)[mask]
        interval_mask = np.isfinite(low) & np.isfinite(high)
        if interval_mask.any():
            coverage = float(np.mean(
                (obs[interval_mask] >= low[interval_mask]) &
                (obs[interval_mask] <= high[interval_mask])
            ))
            interval_width = float(np.mean(high[interval_mask] - low[interval_mask]))
    return {
        "n": n,
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(rss / n)),
        "r2": float(1.0 - rss / denominator) if denominator > 0 else np.nan,
        "calibration_slope": float(slope),
        "calibration_intercept": float(intercept),
        "rss": rss,
        "aic": aic,
        "aicc": aicc,
        "prediction_interval_coverage": coverage,
        "mean_prediction_interval_width": interval_width,
    }


def _closure(frame: pd.DataFrame) -> np.ndarray:
    """Prefer an explicit closure observation, otherwise derive it from widths."""

    for column in OUTPUT_ALIASES["crack_closure_ratio"]:
        if column in frame:
            explicit = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            if np.isfinite(explicit).any():
                return np.clip(explicit, 0.0, 1.0)
    if not {"initial_crack_width_mm", "current_crack_width_mm"}.issubset(frame.columns):
        return np.full(len(frame), np.nan)
    initial = pd.to_numeric(frame["initial_crack_width_mm"], errors="coerce").to_numpy(float)
    current = pd.to_numeric(frame["current_crack_width_mm"], errors="coerce").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        closure = 1.0 - current / initial
    closure[initial <= 0] = np.nan
    return np.clip(closure, 0.0, 1.0)


def _observation(frame: pd.DataFrame, output: str) -> np.ndarray:
    if output == "crack_closure_ratio":
        return _closure(frame)
    for column in OUTPUT_ALIASES[output]:
        if column in frame:
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            if np.isfinite(values).any():
                return values
    return np.full(len(frame), np.nan)


def _predict_outputs(config: ModelConfig, observations: pd.DataFrame) -> pd.DataFrame:
    """Predict every supported response, aligned to input row order."""

    rows = observations.reset_index(drop=True).copy()
    for column in ("initial_crack_width_mm", "wet_hours_per_day", "agent_dosage", "time_d"):
        if column not in rows:
            rows[column] = np.nan
        rows[column] = pd.to_numeric(rows[column], errors="coerce")
    defaults = {
        "initial_crack_width_mm": config.transport.crack_width_mm,
        "wet_hours_per_day": config.environment.wet_hours_per_day,
        "agent_dosage": 1.0,
    }
    group_columns = []
    for column, default in defaults.items():
        key = "_group_{}".format(column)
        rows[key] = rows[column].fillna(default)
        group_columns.append(key)
    result = pd.DataFrame(index=rows.index, columns=OUTPUT_ALIASES.keys(), dtype=float)
    for keys, indexes in rows.groupby(group_columns, sort=False).groups.items():
        width, wet, dosage = (float(value) for value in keys)
        trial = copy.deepcopy(config)
        trial.transport.crack_width_mm = width
        trial.environment.wet_hours_per_day = wet
        trial.kinetics.agent_dosage_multiplier = max(dosage, 0.0)
        subset = rows.loc[indexes]
        finite_times = subset["time_d"].to_numpy(float)
        days = max(float(np.nanmax(finite_times)) if np.isfinite(finite_times).any() else 0.01, 0.01)
        trial.simulation.days = days
        trial.simulation.output_interval_days = max(days / 80.0, 0.01)
        model = simulate_0d(trial).frame
        times = subset["time_d"].fillna(0.0).to_numpy(float)

        def interpolate(column: str) -> np.ndarray:
            return np.interp(times, model["time_d"], model[column])

        cell_volume = interpolate("cell_volume_m3")
        calcite_mol_m3 = interpolate("calcite_mol_m3")
        predictions = {
            "lactate_mM": interpolate("lactate_mol_m3"),
            "calcite_mass_mg": calcite_mol_m3
            * cell_volume * trial.chemistry.calcite_molar_mass_kg_mol * 1.0e6,
            "crack_closure_ratio": interpolate("crack_closure_ratio"),
            "permeability_ratio": interpolate("permeability_ratio"),
            "ph": interpolate("ph"),
            "activation_state": interpolate("activation_state"),
            "cumulative_activity_h": interpolate("cumulative_activity_h"),
        }
        for output, values in predictions.items():
            result.loc[indexes, output] = values
    return result


def _predict_rows(config: ModelConfig, observations: pd.DataFrame) -> np.ndarray:
    """Backward-compatible closure-only prediction helper."""

    return _predict_outputs(config, observations)["crack_closure_ratio"].to_numpy(float)


def _robust_scale(values: np.ndarray, preset: float) -> Tuple[float, str]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) >= 3:
        median = float(np.median(finite))
        mad_scale = 1.4826 * float(np.median(np.abs(finite - median)))
        if mad_scale > 1e-12:
            return mad_scale, "training_mad"
        q25, q75 = np.quantile(finite, [0.25, 0.75])
        iqr_scale = float((q75 - q25) / 1.349)
        if iqr_scale > 1e-12:
            return iqr_scale, "training_iqr"
        std = float(np.std(finite, ddof=1))
        if std > 1e-12:
            return std, "training_std"
    return float(preset), "preset_scale"


def _sigma_for_output(
    frame: pd.DataFrame,
    output: str,
    observed: np.ndarray,
    independent_error: Optional[Mapping[str, float]] = None,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """Apply measurement SD > independent error > robust scale > preset."""

    sigma = np.full(len(frame), np.nan, dtype=float)
    source = np.full(len(frame), "", dtype=object)
    for column in ("{}_sd".format(output), "{}_measurement_sd".format(output)):
        if column in frame:
            candidate = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            mask = np.isfinite(candidate) & (candidate > 0)
            sigma[mask], source[mask] = candidate[mask], "provided:{}.".format(column).rstrip(".")
    if output == "crack_closure_ratio" and "measurement_sd" in frame:
        width_sd = pd.to_numeric(frame["measurement_sd"], errors="coerce").to_numpy(float)
        if "initial_crack_width_mm" in frame:
            initial = pd.to_numeric(frame["initial_crack_width_mm"], errors="coerce").to_numpy(float)
            with np.errstate(divide="ignore", invalid="ignore"):
                candidate = width_sd / initial
        else:
            candidate = width_sd
        mask = np.isnan(sigma) & np.isfinite(candidate) & (candidate > 0)
        sigma[mask], source[mask] = candidate[mask], "provided:measurement_sd_propagated"
    independent = (independent_error or {}).get(output)
    if independent is not None and np.isfinite(independent) and independent > 0:
        mask = np.isnan(sigma)
        sigma[mask], source[mask] = float(independent), "independent_measurement_error"
    robust, robust_source = _robust_scale(observed, PRESET_SCALES[output])
    mask = np.isnan(sigma)
    sigma[mask], source[mask] = robust, robust_source
    return sigma, {
        "minimum": float(np.nanmin(sigma)),
        "median": float(np.nanmedian(sigma)),
        "maximum": float(np.nanmax(sigma)),
        "sources": sorted(set(str(value) for value in source if value)),
    }


def _available_outputs(frame: pd.DataFrame) -> List[str]:
    return [output for output in OUTPUT_ALIASES if np.isfinite(_observation(frame, output)).any()]


def _standardized_residuals(
    config: ModelConfig,
    frame: pd.DataFrame,
    outputs: Sequence[str],
    sigma: Mapping[str, np.ndarray],
) -> Tuple[np.ndarray, Dict[str, float], pd.DataFrame]:
    predictions = _predict_outputs(config, frame)
    residuals, contributions, records = [], {}, []
    for output in outputs:
        observed = _observation(frame, output)
        predicted = predictions[output].to_numpy(float)
        scale = np.asarray(sigma[output], dtype=float)
        mask = np.isfinite(observed) & np.isfinite(predicted) & np.isfinite(scale) & (scale > 0)
        standardized = (predicted[mask] - observed[mask]) / scale[mask]
        residuals.extend(standardized.tolist())
        contributions[output] = float(np.sum(standardized**2))
        for row_index in np.flatnonzero(mask):
            records.append({
                "row_index": int(row_index), "output": output,
                "observed": observed[row_index], "predicted": predicted[row_index],
                "sigma": scale[row_index],
                "standardized_residual": (predicted[row_index] - observed[row_index]) / scale[row_index],
            })
    return np.asarray(residuals, dtype=float), contributions, pd.DataFrame(records)


def _parameter_bounds(config: ModelConfig) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower, upper, initial = [], [], []
    for name in SHARED_PARAMETERS:
        if name in PARAMETER_PROVENANCE:
            bounds = PARAMETER_PROVENANCE[name][2:4]
        else:
            value = _get_parameter(config, name)
            bounds = (max(value / 10.0, 1e-12), value * 10.0)
        lower.append(bounds[0]); upper.append(bounds[1]); initial.append(_get_parameter(config, name))
    return np.asarray(lower), np.asarray(upper), np.asarray(initial)


def _enough_constraints(frame: pd.DataFrame, outputs: Sequence[str], parameter_count: int) -> Tuple[bool, int, int]:
    masks = [np.isfinite(_observation(frame, output)) for output in outputs]
    combined = np.logical_or.reduce(masks) if masks else np.zeros(len(frame), dtype=bool)
    observation_count = int(sum(mask.sum() for mask in masks))
    specimens = int(frame.loc[combined, "specimen_id"].astype(str).nunique()) if combined.any() else 0
    return observation_count >= 2 * parameter_count and specimens >= max(parameter_count, 1), observation_count, specimens


def _fit_empirical_closure(
    frame: pd.DataFrame, sigma: np.ndarray
) -> Tuple[Optional[np.ndarray], Dict[str, object]]:
    observed = _closure(frame)
    times = pd.to_numeric(frame["time_d"], errors="coerce").to_numpy(float)
    mask = np.isfinite(observed) & np.isfinite(times) & np.isfinite(sigma) & (sigma > 0)
    if mask.sum() < 4:
        return None, {"status": "fixed_unavailable", "reason": "fewer than four closure observations"}

    def residual(values: np.ndarray) -> np.ndarray:
        h_inf, rate = values
        prediction = h_inf * (1.0 - np.exp(-rate * times[mask]))
        return (prediction - observed[mask]) / sigma[mask]

    fit = least_squares(residual, np.array([min(max(np.nanmax(observed), 0.05), 1.0), 0.05]),
                        bounds=([0.0, 1e-8], [1.0, 10.0]), max_nfev=200)
    return fit.x, {
        "status": "fitted_on_training_closure",
        "h_inf": float(fit.x[0]), "rate_d_inv": float(fit.x[1]),
        "success": bool(fit.success), "standardized_rss": float(np.sum(fit.fun**2)),
    }


def _fit_stages(
    base: ModelConfig,
    frame: pd.DataFrame,
    outputs: Sequence[str],
    sigma: Mapping[str, np.ndarray],
    initial_log_values: np.ndarray,
    log_lower: np.ndarray,
    log_upper: np.ndarray,
    record_stages: bool = True,
    max_nfev: int = 100,
) -> Tuple[np.ndarray, List[int], List[Dict[str, object]], Dict[str, str]]:
    """Fit each parameter only against its stage-specific observations."""

    current = initial_log_values.copy()
    active_indexes: List[int] = []
    stages, fixed_reasons = [], {}
    for stage_name, candidate_indexes, stage_outputs in STAGES:
        available = [output for output in stage_outputs if output in outputs]
        candidate_indexes = list(candidate_indexes)
        if not candidate_indexes:
            if record_stages:
                stages.append({"stage": stage_name, "status": "evaluated_only",
                               "outputs": available, "parameters": []})
            continue
        enough, observation_count, specimen_count = _enough_constraints(
            frame, available, len(candidate_indexes)
        )
        if not available or not enough:
            reason = "no stage-specific observations" if not available else "insufficient independent constraints"
            for index in candidate_indexes:
                fixed_reasons[SHARED_PARAMETERS[index]] = reason
            if record_stages:
                stages.append({
                    "stage": stage_name, "status": "fixed_to_prior", "outputs": available,
                    "parameters": [SHARED_PARAMETERS[index] for index in candidate_indexes],
                    "observation_count": observation_count, "specimen_count": specimen_count,
                    "reason": reason,
                })
            continue

        def stage_residual(values: np.ndarray) -> np.ndarray:
            trial_values = current.copy()
            trial_values[candidate_indexes] = values
            trial = copy.deepcopy(base)
            for name, value in zip(SHARED_PARAMETERS, np.exp(trial_values)):
                _set_parameter(trial, name, value)
            residuals, _, _ = _standardized_residuals(trial, frame, available, sigma)
            return residuals

        fit = least_squares(
            stage_residual, current[candidate_indexes],
            bounds=(log_lower[candidate_indexes], log_upper[candidate_indexes]),
            max_nfev=max_nfev,
        )
        current[candidate_indexes] = fit.x
        active_indexes.extend(index for index in candidate_indexes if index not in active_indexes)
        jacobian_rank = int(np.linalg.matrix_rank(fit.jac)) if fit.jac.size else 0
        if record_stages:
            stages.append({
                "stage": stage_name, "status": "fitted", "outputs": available,
                "parameters": [SHARED_PARAMETERS[index] for index in candidate_indexes],
                "observation_count": observation_count, "specimen_count": specimen_count,
                "standardized_rss": float(np.sum(fit.fun**2)), "success": bool(fit.success),
                "jacobian_rank": jacobian_rank,
            })
    return current, active_indexes, stages, fixed_reasons


def _empirical_prediction(times: np.ndarray, parameters: Optional[Sequence[float]]) -> np.ndarray:
    if parameters is None:
        return np.full(len(times), np.nan)
    h_inf, rate = parameters
    return float(h_inf) * (1.0 - np.exp(-float(rate) * np.asarray(times, dtype=float)))


def _write_prediction_table(
    path: Path,
    frame: pd.DataFrame,
    config: ModelConfig,
    outputs: Sequence[str],
    sigma: Mapping[str, np.ndarray],
    empirical_parameters: Optional[Sequence[float]],
) -> pd.DataFrame:
    predictions = _predict_outputs(config, frame)
    records = []
    times = pd.to_numeric(frame["time_d"], errors="coerce").to_numpy(float)
    empirical = _empirical_prediction(times, empirical_parameters)
    for output in outputs:
        observed = _observation(frame, output)
        for index in range(len(frame)):
            if not np.isfinite(observed[index]):
                continue
            records.append({
                "row_index": index,
                "dataset_id": frame.iloc[index].get("dataset_id"),
                "specimen_id": frame.iloc[index].get("specimen_id"),
                "split": frame.iloc[index].get("split"),
                "time_d": times[index], "output": output,
                "observed": observed[index], "mechanistic_prediction": predictions.iloc[index][output],
                "sigma": sigma[output][index],
                "standardized_residual": (predictions.iloc[index][output] - observed[index]) / sigma[output][index],
                "empirical_prediction": empirical[index] if output == "crack_closure_ratio" else np.nan,
                "zero_prediction": 0.0 if output == "crack_closure_ratio" else np.nan,
            })
    table = pd.DataFrame(records)
    table.to_csv(path, index=False)
    return table


def _attach_prediction_intervals(
    table: pd.DataFrame,
    observations: pd.DataFrame,
    base: ModelConfig,
    parameter_samples: np.ndarray,
    path: Path,
) -> pd.DataFrame:
    """Propagate bootstrap parameter samples to each observed response."""

    table = table.copy()
    table["prediction_low_95"] = np.nan
    table["prediction_high_95"] = np.nan
    if len(parameter_samples) and not table.empty:
        ensemble = []
        for sample in parameter_samples:
            trial = copy.deepcopy(base)
            for name, value in zip(SHARED_PARAMETERS, sample):
                _set_parameter(trial, name, float(value))
            ensemble.append(_predict_outputs(trial, observations))
        for output in table["output"].dropna().unique():
            values = np.asarray([prediction[output].to_numpy(float) for prediction in ensemble])
            low = np.quantile(values, 0.025, axis=0)
            high = np.quantile(values, 0.975, axis=0)
            selected = table["output"] == output
            source_indexes = table.loc[selected, "row_index"].to_numpy(int)
            table.loc[selected, "prediction_low_95"] = low[source_indexes]
            table.loc[selected, "prediction_high_95"] = high[source_indexes]
    table.to_csv(path, index=False)
    return table


def _write_calibrated_intervals(
    path: Path,
    base: ModelConfig,
    parameter_samples: np.ndarray,
    reference: pd.DataFrame,
    outputs: Sequence[str],
) -> None:
    """Write prospective 7/14/21/28-day parameter-uncertainty intervals."""

    days = np.array([7.0, 14.0, 21.0, 28.0])
    scenarios = pd.DataFrame({"time_d": days})
    for column, default in (
        ("initial_crack_width_mm", base.transport.crack_width_mm),
        ("wet_hours_per_day", base.environment.wet_hours_per_day),
        ("agent_dosage", 1.0),
    ):
        values = pd.to_numeric(reference.get(column, pd.Series(dtype=float)), errors="coerce")
        scenarios[column] = float(values.median()) if values.notna().any() else default
    records = []
    if len(parameter_samples):
        ensembles = {output: [] for output in outputs}
        for sample in parameter_samples:
            trial = copy.deepcopy(base)
            for name, value in zip(SHARED_PARAMETERS, sample):
                _set_parameter(trial, name, float(value))
            prediction = _predict_outputs(trial, scenarios)
            for output in outputs:
                ensembles[output].append(prediction[output].to_numpy(float))
        for output, values in ensembles.items():
            array = np.asarray(values)
            for index, day in enumerate(days):
                records.append({
                    "time_d": day, "output": output,
                    "median": float(np.quantile(array[:, index], 0.5)),
                    "lower_95": float(np.quantile(array[:, index], 0.025)),
                    "upper_95": float(np.quantile(array[:, index], 0.975)),
                    "sample_count": int(len(array)),
                    "evidence_class": "model_prediction",
                })
    pd.DataFrame(records, columns=[
        "time_d", "output", "median", "lower_95", "upper_95",
        "sample_count", "evidence_class",
    ]).to_csv(path, index=False)


def _write_calibration_plots(
    output_dir: Path,
    train_predictions: pd.DataFrame,
    holdout_predictions: pd.DataFrame,
    intervals_path: Path,
) -> None:
    """Create auditable calibration diagnostics without claiming validation."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_dir = output_dir / "calibration_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    outputs = sorted(set(train_predictions.get("output", pd.Series(dtype=str)).dropna()))
    for output in outputs:
        fig, axis = plt.subplots(figsize=(5.4, 4.5))
        values = []
        for table, label, marker in (
            (train_predictions, "Public calibration: train", "o"),
            (holdout_predictions, "Internal specimen holdout", "s"),
        ):
            if table.empty:
                continue
            selected = table.loc[table["output"] == output]
            if selected.empty:
                continue
            axis.scatter(
                selected["observed"], selected["mechanistic_prediction"],
                label=label, marker=marker, alpha=0.8,
            )
            values.extend(selected["observed"].tolist())
            values.extend(selected["mechanistic_prediction"].tolist())
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if len(finite):
            low, high = float(finite.min()), float(finite.max())
            padding = max((high - low) * 0.05, 1e-12)
            axis.plot([low - padding, high + padding], [low - padding, high + padding],
                      color="black", linewidth=1, linestyle="--", label="1:1")
        axis.set_xlabel("Observed ({})".format(output))
        axis.set_ylabel("Mechanistic prediction ({})".format(output))
        axis.set_title("Calibration diagnostic: {}".format(output))
        axis.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(plot_dir / "{}_observed_vs_predicted.png".format(output), dpi=180)
        plt.close(fig)

    if intervals_path.exists():
        intervals = pd.read_csv(intervals_path)
        for output, selected in intervals.groupby("output") if len(intervals) else []:
            fig, axis = plt.subplots(figsize=(5.8, 4.0))
            axis.plot(selected["time_d"], selected["median"], color="#176B87",
                      linewidth=2, label="Bootstrap median")
            axis.fill_between(
                selected["time_d"], selected["lower_95"], selected["upper_95"],
                color="#64CCC5", alpha=0.35, label="95% parameter interval",
            )
            axis.set_xlabel("Time (d)")
            axis.set_ylabel(output)
            axis.set_title("Model prediction: {}".format(output))
            axis.legend(frameon=False, fontsize=8)
            fig.tight_layout()
            fig.savefig(plot_dir / "{}_prediction_interval.png".format(output), dpi=180)
            plt.close(fig)


def calibrate_public(
    data_path: Path,
    output_dir: Path,
    config: Optional[ModelConfig] = None,
    bootstrap_samples: int = 20,
    profile_points: int = 0,
    independent_measurement_error: Optional[Mapping[str, float]] = None,
) -> Dict[str, object]:
    """Calibrate supported outputs in stages and freeze the resulting config."""

    frame = pd.read_csv(data_path)
    if "curation_status" in frame:
        status = frame["curation_status"].fillna("").astype(str)
        if (status == "candidate_only").any():
            raise ValueError(
                "Candidate extraction is not calibration evidence; complete manual curation first"
            )
    required = {"dataset_id", "specimen_id", "split", "time_d"}
    if not required.issubset(frame.columns):
        raise ValueError("Public observation table is missing required columns: {}".format(sorted(required - set(frame))))
    train = frame.loc[frame["split"] == "train"].copy().reset_index(drop=True)
    holdout = frame.loc[frame["split"] == "internal_test"].copy().reset_index(drop=True)
    if train.empty:
        raise ValueError("Public calibration table contains no train rows")
    if set(train["specimen_id"].astype(str)) & set(holdout["specimen_id"].astype(str)):
        raise ValueError("Specimen leakage detected between train and holdout")
    outputs = _available_outputs(train)
    if not outputs:
        raise ValueError("Public calibration table contains no supported observations")

    base = copy.deepcopy(config or ModelConfig())
    lower, upper, initial = _parameter_bounds(base)
    log_lower, log_upper, current = np.log(lower), np.log(upper), np.log(initial)
    train_sigma, holdout_sigma, sigma_audit = {}, {}, {}
    for output in outputs:
        train_sigma[output], audit = _sigma_for_output(
            train, output, _observation(train, output), independent_measurement_error
        )
        holdout_sigma[output], _ = _sigma_for_output(
            holdout, output, _observation(train, output), independent_measurement_error
        ) if len(holdout) else (np.array([]), {})
        sigma_audit[output] = audit

    current, active_indexes, stages, fixed_reasons = _fit_stages(
        base, train, outputs, train_sigma, current, log_lower, log_upper
    )

    if not active_indexes:
        raise ValueError("No parameter stage has sufficient stage-specific observations")

    fitted = copy.deepcopy(base)
    for name, value in zip(SHARED_PARAMETERS, np.exp(current)):
        _set_parameter(fitted, name, value)

    empirical_parameters, empirical_artifact = _fit_empirical_closure(
        train, train_sigma.get("crack_closure_ratio", np.ones(len(train)))
    )
    rng = np.random.RandomState(base.simulation.random_seed)
    bootstrap = []
    specimens = train["specimen_id"].astype(str).unique()
    for _ in range(max(int(bootstrap_samples), 0)):
        sampled_ids = rng.choice(specimens, len(specimens), replace=True)
        blocks = []
        for draw_index, specimen in enumerate(sampled_ids):
            block = train.loc[train["specimen_id"].astype(str) == specimen].copy()
            block["specimen_id"] = "{}__bootstrap_{}".format(specimen, draw_index)
            blocks.append(block)
        sample = pd.concat(blocks, ignore_index=True)
        sample_sigma = {}
        for output in outputs:
            sample_sigma[output], _ = _sigma_for_output(
                sample, output, _observation(train, output), independent_measurement_error
            )
        try:
            values, _, _, _ = _fit_stages(
                base, sample, outputs, sample_sigma, current,
                log_lower, log_upper, record_stages=False, max_nfev=40,
            )
            bootstrap.append(np.exp(values))
        except (RuntimeError, ValueError):
            continue

    output_dir.mkdir(parents=True, exist_ok=True)
    fitted.save(output_dir / "frozen_config.json")
    interval_values = np.asarray(bootstrap) if bootstrap else np.empty((0, len(SHARED_PARAMETERS)))
    boundary_tolerance = 0.01
    parameter_rows, warnings = [], []
    for index, name in enumerate(SHARED_PARAMETERS):
        estimate = _get_parameter(fitted, name)
        log_position = (np.log(estimate) - log_lower[index]) / (log_upper[index] - log_lower[index])
        at_lower, at_upper = log_position <= boundary_tolerance, log_position >= 1.0 - boundary_tolerance
        fitted_flag = index in active_indexes
        ci_low = float(np.quantile(interval_values[:, index], 0.025)) if len(interval_values) else np.nan
        ci_high = float(np.quantile(interval_values[:, index], 0.975)) if len(interval_values) else np.nan
        if fitted_flag and (at_lower or at_upper):
            warnings.append("{} is within 1% of its parameter bound".format(name))
        if fitted_flag and np.isfinite(ci_low) and ci_low > 0 and ci_high / ci_low > 20:
            warnings.append("{} has a bootstrap interval wider than 20-fold".format(name))
        parameter_rows.append({
            "parameter": name, "estimate": estimate, "lower_bound": lower[index], "upper_bound": upper[index],
            "ci_low": ci_low, "ci_high": ci_high, "fitted": fitted_flag,
            "fixed_reason": fixed_reasons.get(name, ""), "at_lower_boundary": bool(at_lower),
            "at_upper_boundary": bool(at_upper),
        })
    pd.DataFrame(parameter_rows).to_csv(output_dir / "fitted_parameters.csv", index=False)

    if len(interval_values):
        bootstrap_frame = pd.DataFrame(interval_values, columns=SHARED_PARAMETERS)
        bootstrap_frame.to_csv(output_dir / "bootstrap_parameters.csv", index=False)
        correlation = bootstrap_frame[[SHARED_PARAMETERS[index] for index in active_indexes]].corr()
        correlation.to_csv(output_dir / "parameter_correlation.csv")
        for row_index in range(len(correlation)):
            for column_index in range(row_index + 1, len(correlation)):
                if abs(correlation.iloc[row_index, column_index]) >= 0.95:
                    warnings.append("high bootstrap correlation between {} and {}".format(
                        correlation.index[row_index], correlation.columns[column_index]
                    ))
    else:
        pd.DataFrame(index=[SHARED_PARAMETERS[index] for index in active_indexes]).to_csv(
            output_dir / "parameter_correlation.csv"
        )
    _write_calibrated_intervals(
        output_dir / "calibrated_prediction_intervals.csv",
        fitted, interval_values, train, outputs,
    )

    profile_rows = []
    if profile_points >= 3:
        for fixed_index in active_indexes:
            profile_outputs = [
                output for output in outputs
                if fixed_index in OUTPUT_PARAMETER_INDEXES.get(output, ())
            ]
            stage_indexes = sorted({
                index for output in profile_outputs
                for index in OUTPUT_PARAMETER_INDEXES.get(output, ())
                if index in active_indexes
            })
            free_indexes = [index for index in stage_indexes if index != fixed_index]
            for fixed_value in np.linspace(log_lower[fixed_index], log_upper[fixed_index], profile_points):
                def profile_residual(free_values: np.ndarray) -> np.ndarray:
                    values = current.copy()
                    values[fixed_index] = fixed_value
                    values[free_indexes] = free_values
                    trial = copy.deepcopy(base)
                    for name, value in zip(SHARED_PARAMETERS, np.exp(values)):
                        _set_parameter(trial, name, value)
                    residuals, _, _ = _standardized_residuals(
                        trial, train, profile_outputs, train_sigma
                    )
                    return residuals
                if free_indexes:
                    profile_fit = least_squares(
                        profile_residual, current[free_indexes],
                        bounds=(log_lower[free_indexes], log_upper[free_indexes]), max_nfev=50,
                    )
                    cost, success = float(np.sum(profile_fit.fun**2)), bool(profile_fit.success)
                else:
                    values = profile_residual(np.array([]))
                    cost, success = float(np.sum(values**2)), True
                profile_rows.append({
                    "parameter": SHARED_PARAMETERS[fixed_index], "fixed_value": float(np.exp(fixed_value)),
                    "standardized_rss": cost, "delta_rss": np.nan, "success": success,
                })
        profile_frame = pd.DataFrame(profile_rows)
        profile_frame["delta_rss"] = profile_frame.groupby("parameter")["standardized_rss"].transform(
            lambda values: values - values.min()
        )
        profile_frame.to_csv(output_dir / "profile_likelihood.csv", index=False)
    else:
        pd.DataFrame(columns=["parameter", "fixed_value", "standardized_rss", "delta_rss", "success"]).to_csv(
            output_dir / "profile_likelihood.csv", index=False
        )

    train_table = _write_prediction_table(
        output_dir / "train_predictions.csv", train, fitted, outputs, train_sigma, empirical_parameters
    )
    holdout_table = _write_prediction_table(
        output_dir / "internal_test_predictions.csv", holdout, fitted, outputs, holdout_sigma, empirical_parameters
    ) if len(holdout) else pd.DataFrame()
    train_table = _attach_prediction_intervals(
        train_table, train, fitted, interval_values,
        output_dir / "train_predictions.csv",
    )
    if len(holdout):
        holdout_table = _attach_prediction_intervals(
            holdout_table, holdout, fitted, interval_values,
            output_dir / "internal_test_predictions.csv",
        )
    _write_calibration_plots(
        output_dir, train_table, holdout_table,
        output_dir / "calibrated_prediction_intervals.csv",
    )
    _, contributions, _ = _standardized_residuals(fitted, train, outputs, train_sigma)
    total_contribution = max(sum(contributions.values()), 1e-30)
    contribution_report = {
        output: {"standardized_rss": value, "fraction": value / total_contribution}
        for output, value in contributions.items()
    }
    active_count = len(active_indexes)
    metrics = {"train": {}, "internal_test": {}}
    for split_name, table in (("train", train_table), ("internal_test", holdout_table)):
        for output in outputs:
            output_parameter_count = len(
                set(OUTPUT_PARAMETER_INDEXES.get(output, ())) & set(active_indexes)
            )
            selected = table.loc[table["output"] == output] if len(table) else pd.DataFrame()
            if selected.empty:
                metrics[split_name][output] = _metrics(
                    np.array([]), np.array([]), output_parameter_count
                )
            else:
                metrics[split_name][output] = _metrics(
                    selected["observed"], selected["mechanistic_prediction"], output_parameter_count,
                    selected["prediction_low_95"], selected["prediction_high_95"],
                )
    closure_train = train_table.loc[train_table["output"] == "crack_closure_ratio"]
    closure_test = holdout_table.loc[holdout_table["output"] == "crack_closure_ratio"] if len(holdout_table) else pd.DataFrame()
    baseline_metrics = {"empirical_parameters": empirical_artifact, "train": {}, "internal_test": {}}
    for split_name, table in (("train", closure_train), ("internal_test", closure_test)):
        if len(table):
            baseline_metrics[split_name] = {
                "mechanistic": _metrics(table["observed"], table["mechanistic_prediction"], 1),
                "empirical": _metrics(table["observed"], table["empirical_prediction"], 2),
                "zero": _metrics(table["observed"], table["zero_prediction"], 0),
            }

    stage_ranks = {
        stage["stage"]: {
            "jacobian_rank": stage.get("jacobian_rank", 0),
            "parameter_count": len(stage.get("parameters", [])),
        }
        for stage in stages if stage.get("status") == "fitted"
    }
    jacobian_rank = sum(item["jacobian_rank"] for item in stage_ranks.values())
    for stage_name, rank in stage_ranks.items():
        if rank["jacobian_rank"] < rank["parameter_count"]:
            warnings.append("{} calibration Jacobian is rank deficient ({}/{})".format(
                stage_name, rank["jacobian_rank"], rank["parameter_count"]
            ))
    joint_residuals, _, _ = _standardized_residuals(fitted, train, outputs, train_sigma)
    standardized_rss = float(np.sum(joint_residuals**2))
    standardized_n = int(len(joint_residuals))
    joint_aic = float(standardized_n * np.log(max(standardized_rss / standardized_n, 1e-30)) + 2 * active_count)
    joint_aicc = (
        float(joint_aic + 2 * active_count * (active_count + 1) / (standardized_n - active_count - 1))
        if standardized_n > active_count + 1 else np.nan
    )
    calibration_metrics = {
        "evidence_class": "public_calibration_data",
        "outputs": outputs,
        "sigma_audit": sigma_audit,
        "per_output": metrics,
        "joint_objective": {
            "standardized_residual_count": standardized_n,
            "standardized_rss": standardized_rss,
            "parameter_count": active_count,
            "aic": joint_aic, "aicc": joint_aicc,
            "loss_contributions": contribution_report,
        },
        "baseline_comparison": baseline_metrics,
    }
    (output_dir / "calibration_metrics.json").write_text(
        json.dumps(calibration_metrics, indent=2), encoding="utf-8"
    )

    artifact = {
        "evidence_class": "public_calibration_data",
        "project_experiment_rows": 0,
        "dataset_ids": sorted(frame["dataset_id"].dropna().unique().tolist()),
        "frozen_config_sha256": _digest_config(fitted),
        "shared_parameters": list(SHARED_PARAMETERS),
        "fitted_parameters": [SHARED_PARAMETERS[index] for index in active_indexes],
        "fixed_to_prior": fixed_reasons,
        "calibration_stages": stages,
        "empirical_baseline": empirical_artifact,
        "bootstrap_successes": int(len(interval_values)),
        "profile_likelihood_points": int(profile_points),
        "identifiability": {"jacobian_rank": jacobian_rank, "active_parameter_count": active_count,
                            "stage_ranks": stage_ranks, "warnings": sorted(set(warnings))},
        "joint_objective": calibration_metrics["joint_objective"],
    }
    (output_dir / "frozen_run.json").write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    return artifact


def validate_external(data_path: Path, frozen_run_dir: Path, output_dir: Path) -> Dict[str, object]:
    """Evaluate a frozen model. This function intentionally contains no optimizer."""

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
        prediction_low = np.quantile(ensemble_values, 0.025, axis=0)
        prediction_high = np.quantile(ensemble_values, 0.975, axis=0)
    else:
        prediction_low = np.full(len(frame), np.nan)
        prediction_high = np.full(len(frame), np.nan)
    times = pd.to_numeric(frame["time_d"], errors="coerce").to_numpy(float)
    empirical_artifact = artifact.get("empirical_baseline", {})
    empirical_parameters = None
    if empirical_artifact.get("status") == "fitted_on_training_closure":
        empirical_parameters = [empirical_artifact["h_inf"], empirical_artifact["rate_d_inv"]]
    empirical = _empirical_prediction(times, empirical_parameters)
    zero = np.zeros(len(frame))
    active_count = len(artifact.get("fitted_parameters", []))
    results = {
        "evidence_class": "external_validation_data",
        "parameters_refitted": False,
        "frozen_config_sha256": artifact["frozen_config_sha256"],
        "mechanistic": _metrics(observed, mechanism, active_count, prediction_low, prediction_high),
        "first_order_baseline": _metrics(observed, empirical, 2),
        "zero_mineralization_baseline": _metrics(observed, zero, 0),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "observed_closure": observed, "mechanistic_prediction": mechanism,
        "prediction_low_95": prediction_low, "prediction_high_95": prediction_high,
        "first_order_baseline": empirical, "zero_baseline": zero,
    }).to_csv(output_dir / "predictions.csv", index=False)
    (output_dir / "external_validation.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return results


def fit_measurement_error(data_path: Path, output_dir: Path) -> Dict[str, object]:
    frame = pd.read_csv(data_path)
    numeric = frame.select_dtypes(include=[np.number])
    reference = next((name for name in numeric if "reference" in name.lower() or "manual" in name.lower()), None)
    prediction = next((name for name in numeric if name != reference and
                       ("predict" in name.lower() or "model" in name.lower())), None)
    if reference is None or prediction is None:
        raise ValueError("Measurement table needs identifiable reference/manual and prediction columns")
    residual = numeric[prediction] - numeric[reference]
    result = {
        "evidence_class": "measurement_error_data",
        "rows": int(residual.notna().sum()), "bias": float(residual.mean()),
        "measurement_sd": float(residual.std(ddof=1)), "kinetic_parameters_fitted": False,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "measurement_error.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result

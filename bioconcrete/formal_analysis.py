"""Standards-based Morris and Sobol analysis using SALib."""

from __future__ import annotations

import copy
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from .analysis import _set_parameter
from .config import ModelConfig
from .model import simulate_0d


FORMAL_PARAMETERS = (
    "kinetics.activity_multiplier",
    "kinetics.capsule_release_s",
    "kinetics.response_delay_h",
    "kinetics.basal_leak_fraction",
    "kinetics.k_oxygen_mol_m3",
    "chemistry.wall_deposition_fraction",
    "chemistry.calcite_rate_mol_m3_s",
    "transport.crack_width_mm",
    "environment.wet_hours_per_day",
)

FORMAL_BOUNDS = (
    (0.5, 5.0), (1e-7, 3e-5), (0.0, 24.0), (0.0, 0.10),
    (0.005, 0.20), (0.05, 1.0), (1e-7, 1e-2), (0.1, 0.5), (6.0, 24.0),
)


def _problem(parameters: Sequence[str] = FORMAL_PARAMETERS) -> Dict[str, object]:
    indexes = [FORMAL_PARAMETERS.index(name) for name in parameters]
    return {"num_vars": len(parameters), "names": list(parameters),
            "bounds": [list(FORMAL_BOUNDS[index]) for index in indexes]}


def _evaluate_matrix(matrix: np.ndarray, problem: Dict[str, object], config: ModelConfig) -> np.ndarray:
    output = np.empty(len(matrix))
    for row_index, row in enumerate(matrix):
        trial = copy.deepcopy(config)
        trial.simulation.output_interval_days = trial.simulation.days
        for name, value in zip(problem["names"], row):
            _set_parameter(trial, name, float(value))
        output[row_index] = simulate_0d(trial).summary["mean_crack_closure_ratio"]
    return output


def _config_from_dict(raw: Dict[str, object]) -> ModelConfig:
    config = ModelConfig()
    for section, values in raw.items():
        for name, value in values.items():
            setattr(getattr(config, section), name, value)
    config.validate()
    return config


def _sample_id(method: str, index: int, row: np.ndarray) -> str:
    payload = method + "|" + str(index) + "|" + ",".join("{:.17g}".format(value) for value in row)
    return hashlib.sha256(payload.encode()).hexdigest()[:20]


def _evaluate_sample(payload: Tuple[str, int, np.ndarray, Dict[str, object], Dict[str, object]]) -> Dict[str, object]:
    method, index, row, problem, raw = payload
    config = _config_from_dict(raw)
    config.simulation.output_interval_days = config.simulation.days
    for name, value in zip(problem["names"], row):
        _set_parameter(config, name, float(value))
    value = simulate_0d(config).summary["mean_crack_closure_ratio"]
    return {"sample_id": _sample_id(method, index, row), "method": method,
            "sample_index": index, **{name: float(value) for name, value in zip(problem["names"], row)},
            "response": float(value)}


def evaluate_resumable(
    matrices: Dict[str, np.ndarray], problem: Dict[str, object], config: ModelConfig,
    output_dir: Path, workers: int = 1, resume: bool = False,
) -> Tuple[Dict[str, np.ndarray], pd.DataFrame]:
    """Evaluate deterministic sensitivity designs with resume and failure traces."""

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "sensitivity_samples.csv"
    existing = pd.read_csv(path) if resume and path.exists() else pd.DataFrame()
    rows = existing.to_dict("records") if len(existing) else []
    completed = set(existing.get("sample_id", pd.Series(dtype=str)).astype(str))
    tasks = []
    raw = config.to_dict()
    for method, matrix in matrices.items():
        for index, row in enumerate(matrix):
            if _sample_id(method, index, row) not in completed:
                tasks.append((method, index, row, problem, raw))
    failures = []

    pending_writes = {"count": 0}

    def flush() -> None:
        frame = pd.DataFrame(rows).sort_values(["method", "sample_index"])
        temporary = path.with_suffix(".tmp")
        frame.to_csv(temporary, index=False)
        temporary.replace(path)
        pending_writes["count"] = 0

    def record(row: Dict[str, object]) -> None:
        rows.append(row)
        pending_writes["count"] += 1
        if pending_writes["count"] >= 16:
            flush()

    if workers <= 1:
        for task in tasks:
            try:
                record(_evaluate_sample(task))
            except Exception as error:
                failures.append({"sample_id": _sample_id(task[0], task[1], task[2]), "error": str(error)})
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_evaluate_sample, task): task for task in tasks}
            for future in as_completed(futures):
                try:
                    record(future.result())
                except Exception as error:
                    task = futures[future]
                    failures.append({"sample_id": _sample_id(task[0], task[1], task[2]), "error": str(error)})
    if pending_writes["count"] or not path.exists():
        flush()
    frame = pd.DataFrame(rows).drop_duplicates("sample_id", keep="last")
    temporary = path.with_suffix(".tmp")
    frame.sort_values(["method", "sample_index"]).to_csv(temporary, index=False)
    temporary.replace(path)
    pd.DataFrame(failures, columns=["sample_id", "error"]).to_csv(output_dir / "failed_samples.csv", index=False)
    outputs = {}
    for method, matrix in matrices.items():
        selected = frame.loc[frame["method"] == method].sort_values("sample_index")
        if len(selected) != len(matrix):
            raise RuntimeError("{} sensitivity design incomplete: {}/{}".format(method, len(selected), len(matrix)))
        outputs[method] = selected["response"].to_numpy(float)
    return outputs, frame


def formal_sensitivity(output_dir: Path, config: Optional[ModelConfig] = None,
                       samples: int = 256, seed: int = 2026, workers: int = 1,
                       resume: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run independent-trajectory Morris and Saltelli/Sobol designs."""

    try:
        from SALib.analyze import morris as morris_analyze
        from SALib.analyze import sobol as sobol_analyze
        from SALib.sample import morris as morris_sample
        from SALib.sample import saltelli
    except ImportError as error:
        raise RuntimeError("formal-sensitivity requires SALib 1.4.x; install the project dependencies") from error
    base = copy.deepcopy(config or ModelConfig())
    problem = _problem()
    samples = max(int(samples), 16)
    np.random.seed(seed)
    morris_x = morris_sample.sample(problem, N=max(16, samples // 4), num_levels=4,
                                    optimal_trajectories=None, seed=seed)
    sobol_x = saltelli.sample(problem, samples, calc_second_order=False)
    evaluated, sample_frame = evaluate_resumable(
        {"morris": morris_x, "sobol": sobol_x}, problem, base, output_dir, workers, resume
    )
    morris_y = evaluated["morris"]
    morris_result = morris_analyze.analyze(problem, morris_x, morris_y, num_levels=4,
                                           num_resamples=500, conf_level=.95,
                                           print_to_console=False, seed=seed)
    morris = pd.DataFrame({
        "parameter": problem["names"], "mu": morris_result["mu"],
        "mu_star": morris_result["mu_star"], "mu_star_conf": morris_result["mu_star_conf"],
        "sigma": morris_result["sigma"],
    })
    sobol_y = evaluated["sobol"]
    sobol_result = sobol_analyze.analyze(problem, sobol_y, calc_second_order=False,
                                         num_resamples=500, conf_level=.95,
                                         print_to_console=False, seed=seed)
    sobol = pd.DataFrame({
        "parameter": problem["names"], "S1": sobol_result["S1"],
        "S1_conf": sobol_result["S1_conf"], "ST": sobol_result["ST"],
        "ST_conf": sobol_result["ST_conf"],
    })
    sobol["converged"] = ~(
        (sobol["S1"] < -sobol["S1_conf"]) | (sobol["S1"] > 1 + sobol["S1_conf"]) |
        (sobol["ST"] < -sobol["ST_conf"]) | (sobol["ST"] > 1 + sobol["ST_conf"])
    )
    block_size = problem["num_vars"] + 2
    convergence_rows = []
    levels = sorted(set(max(16, samples // divisor) for divisor in (4, 2, 1)))
    for level in levels:
        level_result = sobol_analyze.analyze(
            problem, sobol_y[:level * block_size], calc_second_order=False,
            num_resamples=200, conf_level=.95, print_to_console=False, seed=seed,
        )
        for index, name in enumerate(problem["names"]):
            convergence_rows.append({
                "base_samples": level, "parameter": name,
                "S1": level_result["S1"][index], "ST": level_result["ST"][index],
                "S1_conf": level_result["S1_conf"][index], "ST_conf": level_result["ST_conf"][index],
            })
    output_dir.mkdir(parents=True, exist_ok=True)
    morris.to_csv(output_dir / "morris_indices.csv", index=False)
    sobol.to_csv(output_dir / "sobol_indices.csv", index=False)
    pd.DataFrame(convergence_rows).to_csv(output_dir / "sobol_convergence.csv", index=False)
    metadata = {"implementation": "SALib", "base_samples": samples,
                "sobol_evaluations": int(len(sobol_x)), "morris_evaluations": int(len(morris_x)),
                "indices_clipped": False, "workers": workers, "resume": resume,
                "completed_evaluations": len(sample_frame),
                "warning": "Rows marked unconverged must not be interpreted."}
    (output_dir / "formal_sensitivity.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return morris, sobol

"""Standards-based Morris and Sobol analysis using SALib."""

from __future__ import annotations

import copy
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


def formal_sensitivity(output_dir: Path, config: Optional[ModelConfig] = None,
                       samples: int = 256, seed: int = 2026) -> Tuple[pd.DataFrame, pd.DataFrame]:
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
    morris_y = _evaluate_matrix(morris_x, problem, base)
    morris_result = morris_analyze.analyze(problem, morris_x, morris_y, num_levels=4,
                                           num_resamples=500, conf_level=.95,
                                           print_to_console=False, seed=seed)
    morris = pd.DataFrame({
        "parameter": problem["names"], "mu": morris_result["mu"],
        "mu_star": morris_result["mu_star"], "mu_star_conf": morris_result["mu_star_conf"],
        "sigma": morris_result["sigma"],
    })
    sobol_x = saltelli.sample(problem, samples, calc_second_order=False)
    sobol_y = _evaluate_matrix(sobol_x, problem, base)
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
                "indices_clipped": False, "warning": "Rows marked unconverged must not be interpreted."}
    (output_dir / "formal_sensitivity.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return morris, sobol

"""Numerical and physical acceptance checks."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict

from .config import ModelConfig
from .model import simulate_0d, simulate_1d, simulate_2d


def _relative_change(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), abs(second), 1e-12)


def run_validation(output_dir: Path, config: ModelConfig, full: bool = False) -> Dict[str, object]:
    """Run conservation, limiting-case, and discretization checks."""

    output_dir.mkdir(parents=True, exist_ok=True)
    checks: Dict[str, object] = {}

    closed = copy.deepcopy(config)
    closed.simulation.closed_system = True
    closed.simulation.days = config.simulation.days if full else min(config.simulation.days, 1.0)
    closed_result = simulate_0d(closed)
    balance = closed_result.diagnostics["relative_balance_change"]
    checks["closed_carbon_balance"] = {"value": balance["carbon"], "limit": 0.005, "passed": balance["carbon"] < 0.005}
    checks["closed_calcium_balance"] = {"value": balance["calcium"], "limit": 0.005, "passed": balance["calcium"] < 0.005}
    checks["ammonia_free"] = {"value": closed_result.summary["ammonium_mol_m3_max"], "limit": 0.0, "passed": closed_result.summary["ammonium_mol_m3_max"] == 0.0}
    checks["nonnegative"] = {"passed": bool(closed_result.diagnostics["nonnegative"])}

    no_source = copy.deepcopy(closed)
    no_source.kinetics.capsule_calcium_lactate_mol_m3 = 0.0
    no_source.kinetics.spore_density_rel = 0.0
    no_source.kinetics.active_density_rel = 0.0
    no_source_result = simulate_0d(no_source)
    checks["no_source_no_mineralization"] = {
        "value": no_source_result.summary["calcite_mol_m3_mean"],
        "limit": 1e-10,
        "passed": no_source_result.summary["calcite_mol_m3_mean"] < 1e-10,
    }

    coarse_time = copy.deepcopy(config)
    fine_time = copy.deepcopy(config)
    if not full:
        coarse_time.simulation.days = fine_time.simulation.days = min(config.simulation.days, 1.0)
    fine_time.simulation.reaction_step_h = coarse_time.simulation.reaction_step_h / 2.0
    coarse_result, fine_result = simulate_0d(coarse_time), simulate_0d(fine_time)
    time_change = _relative_change(coarse_result.summary["mean_healing_ratio"], fine_result.summary["mean_healing_ratio"])
    checks["time_step_convergence"] = {"value": time_change, "limit": 0.05, "passed": time_change < 0.05}

    spatial_days = config.simulation.days if full else min(config.simulation.days, 0.1)
    one_coarse, one_fine = copy.deepcopy(config), copy.deepcopy(config)
    one_coarse.simulation.days = one_fine.simulation.days = spatial_days
    if not full:
        one_coarse.transport.nx_1d, one_fine.transport.nx_1d = 9, 17
    else:
        one_fine.transport.nx_1d = 2 * one_coarse.transport.nx_1d - 1
    one_a, one_b = simulate_1d(one_coarse), simulate_1d(one_fine)
    one_change = _relative_change(one_a.summary["mean_healing_ratio"], one_b.summary["mean_healing_ratio"])
    checks["one_dimensional_grid_convergence"] = {"value": one_change, "limit": 0.05, "passed": one_change < 0.05}

    two_coarse, two_fine = copy.deepcopy(config), copy.deepcopy(config)
    two_coarse.simulation.days = two_fine.simulation.days = spatial_days
    if not full:
        two_coarse.transport.nx_2d, two_coarse.transport.ny_2d = 7, 3
        two_fine.transport.nx_2d, two_fine.transport.ny_2d = 13, 5
    else:
        two_fine.transport.nx_2d = 2 * two_coarse.transport.nx_2d - 1
        two_fine.transport.ny_2d = 2 * two_coarse.transport.ny_2d - 1
    two_a, two_b = simulate_2d(two_coarse), simulate_2d(two_fine)
    two_change = _relative_change(two_a.summary["mean_healing_ratio"], two_b.summary["mean_healing_ratio"])
    checks["two_dimensional_grid_convergence"] = {"value": two_change, "limit": 0.05, "passed": two_change < 0.05}

    report = {
        "mode": "full" if full else "quick_screen",
        "all_passed": all(item.get("passed", False) for item in checks.values()),
        "checks": checks,
        "interpretation": "A quick spatial screen does not replace full-duration experimental validation.",
    }
    (output_dir / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report

"""Gate-D verification for coupled three-dimensional simulations."""
from __future__ import annotations

from copy import deepcopy
import csv
import json
import math
import os
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

from .config import ModelConfig
from .grid_3d import rectangular_grid_3d
from .io_3d import config_hash
from .model import _initial_state, _reaction_step
from .model_3d import EVIDENCE_LABEL, simulate_3d
from .reaction_kernel import reaction_step_cells
from .transport_3d import transport_step_3d


METRICS = (
    "calcite_mol", "area_weighted_closure", "maximum_local_closure",
    "open_volume_m3", "oxygen_penetration_depth_m",
)
NOISE_FLOORS = {
    "calcite_mol": 1e-15,
    "area_weighted_closure": 1e-10,
    "maximum_local_closure": 1e-10,
    "open_volume_m3": 1e-15,
    "oxygen_penetration_depth_m": 1e-8,
}


def validation_config(config: ModelConfig = None) -> ModelConfig:
    """Freeze the nonzero one-day numerical verification scenario."""
    c = deepcopy(config or ModelConfig())
    c.simulation.days = 1.0
    c.simulation.output_interval_days = 1.0
    c.output_3d.save_every_days = 1.0
    c.simulation.reaction_step_h = 6.0
    c.simulation.closed_system = True
    c.simulation.random_seed = 2026
    c.environment.exposure = "continuous"
    c.environment.oxygen_initial_mol_m3 = c.environment.oxygen_boundary_mol_m3
    c.environment.inorganic_carbon_boundary_mol_m3 = 0.15
    c.environment.inorganic_carbon_initial_mol_m3 = 0.15
    c.kinetics.gate_logic = "static_suitability"
    c.kinetics.response_delay_h = 24.0
    c.kinetics.activation_duration_h = 24.0
    c.kinetics.signal_relaxation_h = 24.0
    c.kinetics.basal_leak_fraction = 1.0
    c.kinetics.capsule_calcium_lactate_mol_m3 = 80.0
    c.kinetics.capsule_release_s = 1.0e-7
    c.kinetics.csh_release_s = 1.0e-7
    c.kinetics.maximum_growth_s = 1.0e-7
    c.kinetics.effective_kcat_s = 1.0e-5
    c.chemistry.portlandite_mol_m3 = 50.0
    c.chemistry.calcite_rate_mol_m3_s = 1.0e-7
    c.chemistry.portlandite_dissolution_s = 1.0e-8
    c.geometry_3d.capsule_depth_mode = "uniform"
    # Grid convergence must sample the same resolved continuous inventory field,
    # not three differently under-resolved narrow capsule kernels.
    c.geometry_3d.capsule_count = 8
    c.geometry_3d.capsule_spread_x_mm = c.transport.crack_length_mm
    c.geometry_3d.capsule_spread_y_mm = c.transport.crack_width_mm
    c.geometry_3d.capsule_spread_z_mm = c.transport.crack_depth_mm
    c.solver_3d.splitting_scheme = "strang"
    c.solver_3d.linear_solver = "auto"
    # Moderate sparse blocks amortize the few exceptionally stiff cells without
    # constructing one global system. SciPy BDF is kept serial to avoid BLAS
    # oversubscription; the parallel path has a separate equivalence test.
    c.solver_3d.reaction_batch_cells = 16
    c.solver_3d.reaction_workers = max(1, min(8, os.cpu_count() or 1))
    c.solver_3d.reaction_parallel_backend = "process"
    c.validate()
    return c


def _relative_error(first: float, second: float, floor: float) -> Dict[str, object]:
    scale = max(abs(first), abs(second))
    if scale < floor:
        return {"first": first, "second": second, "absolute_error": abs(first-second),
                "relative_error": None, "excited": False, "passed": True}
    error = abs(first-second) / scale
    return {"first": first, "second": second, "absolute_error": abs(first-second),
            "relative_error": error, "excited": True, "passed": error < 0.05}


def dimensional_reduction_checks(config: ModelConfig = None) -> Dict[str, object]:
    c = validation_config(config)
    c.simulation.days = 0.01
    c.simulation.reaction_step_h = 0.24
    dt = c.simulation.days * 86400.0
    local = _initial_state(c, np.ones(1), "0d")
    old = _reaction_step(local.copy(), 0.0, dt, c, None)
    shared = reaction_step_cells(local.copy(), 0.0, dt, c, batch_size=1)
    zero_d_error = float(np.max(np.abs(old-shared) / np.maximum(np.abs(old), 1e-12)))

    c.geometry_3d.nx, c.geometry_3d.ny, c.geometry_3d.nz = 7, 3, 5
    grid = rectangular_grid_3d(c)
    x_profile = np.linspace(0.0, 1.0, grid.shape[2])
    field = np.broadcast_to(x_profile, grid.shape).copy()
    transported, _ = transport_step_3d(field, grid, 10.0, 1e-9)
    yz_uniform_error = float(np.max(np.ptp(transported, axis=(0, 1))))

    xy = np.arange(grid.shape[1] * grid.shape[2], dtype=float).reshape(grid.shape[1], grid.shape[2])
    z_uniform = np.broadcast_to(xy, grid.shape).copy()
    z_result, _ = transport_step_3d(z_uniform, grid, 10.0, 1e-9)
    z_uniform_error = float(np.max(np.ptp(z_result, axis=0)))

    xz = np.arange(grid.shape[0] * grid.shape[2], dtype=float).reshape(grid.shape[0], grid.shape[2])
    y_uniform = np.broadcast_to(xz[:, None, :], grid.shape).copy()
    y_result, _ = transport_step_3d(y_uniform, grid, 10.0, 1e-9)
    y_uniform_error = float(np.max(np.ptp(y_result, axis=1)))
    return {
        "single_voxel_to_0d": {"relative_error": zero_d_error, "limit": 0.01,
                                "passed": zero_d_error < 0.01},
        "yz_uniform_to_1d_transport": {"absolute_error": yz_uniform_error,
                                        "limit": 1e-6, "passed": yz_uniform_error < 1e-6},
        "z_uniform_to_legacy_xy_2d_transport": {"absolute_error": z_uniform_error,
                                                 "limit": 1e-6, "passed": z_uniform_error < 1e-6},
        "y_uniform_to_xz_2p5d_transport": {"absolute_error": y_uniform_error,
                                           "limit": 1e-6, "passed": y_uniform_error < 1e-6},
        "note": "Legacy 2D is (x,y); production 2.5D is (x,z). The adapters are intentionally distinct.",
    }


def _run_grid(config: ModelConfig, grid: Tuple[int, int, int]) -> Dict[str, float]:
    c = deepcopy(config)
    c.geometry_3d.nx, c.geometry_3d.ny, c.geometry_3d.nz = grid
    return simulate_3d(c).summary


def _cached_run(output: Path, name: str, config: ModelConfig,
                grid: Tuple[int, int, int]) -> Dict[str, float]:
    cache = output / "scenarios" / (name + ".json")
    if cache.exists():
        payload = json.loads(cache.read_text(encoding="utf-8"))
        probe = deepcopy(config)
        probe.geometry_3d.nx, probe.geometry_3d.ny, probe.geometry_3d.nz = grid
        if payload.get("config_hash") == config_hash(probe):
            return payload["summary"]
    c = deepcopy(config)
    c.geometry_3d.nx, c.geometry_3d.ny, c.geometry_3d.nz = grid
    result = simulate_3d(c)
    cache.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache.with_suffix(".json.tmp")
    temporary.write_text(json.dumps({"config_hash": config_hash(c), "summary": result.summary,
                                     "performance": result.performance}, indent=2), encoding="utf-8")
    temporary.replace(cache)
    return result.summary


def _comparison(first: Dict[str, float], second: Dict[str, float]) -> Dict[str, object]:
    return {name: _relative_error(first[name], second[name], NOISE_FLOORS[name]) for name in METRICS}


def run_validation_3d(output: Path, config: ModelConfig = None, full: bool = False) -> Dict[str, object]:
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    base = validation_config(config)
    reduction = dimensional_reduction_checks(base)
    smoke = deepcopy(base)
    smoke.simulation.days = 0.01
    smoke.output_3d.save_every_days = 0.01
    smoke.simulation.reaction_step_h = 0.24
    smoke.geometry_3d.nx, smoke.geometry_3d.ny, smoke.geometry_3d.nz = 9, 3, 5
    smoke_result = simulate_3d(smoke)
    conservation = smoke_result.diagnostics["conservation"]
    conservation_passed = all(conservation[key]["relative_error"] < 0.005
                              for key in ("carbon_mol", "calcium_mol"))
    reduction_passed = all(value.get("passed", True) for value in reduction.values()
                           if isinstance(value, dict))

    grid_results = time_results = {}
    grid_comparison = time_comparison = {}
    if full:
        grids = {"coarse": (21, 3, 9), "medium": (41, 5, 17), "fine": (81, 7, 33)}
        grid_base = deepcopy(base)
        grid_base.simulation.reaction_step_h = 24.0
        grid_results = {name: _cached_run(output, "grid_" + name, grid_base, shape)
                        for name, shape in grids.items()}
        grid_comparison = _comparison(grid_results["medium"], grid_results["fine"])
        for label, hours in (("dt", 6.0), ("dt_half", 3.0), ("dt_quarter", 1.5)):
            c = deepcopy(base); c.simulation.reaction_step_h = hours
            time_results[label] = _cached_run(output, "time_" + label, c, (31, 3, 11))
        time_comparison = _comparison(time_results["dt_half"], time_results["dt_quarter"])
    convergence_passed = bool(full and
        all(item["passed"] for item in grid_comparison.values()) and
        all(item["passed"] for item in time_comparison.values()))
    gate_passed = bool(full and conservation_passed and reduction_passed and
                       convergence_passed and smoke_result.diagnostics["finite"] and
                       smoke_result.diagnostics["nonnegative"])
    report = {
        "model_version": "v0.6.0-development", "evidence_label": EVIDENCE_LABEL,
        "full": full, "gate_d_passed": gate_passed,
        "checks": {"conservation_passed": conservation_passed,
                   "reduction_passed": reduction_passed,
                   "convergence_passed": convergence_passed,
                   "finite": smoke_result.diagnostics["finite"],
                   "nonnegative": smoke_result.diagnostics["nonnegative"]},
        "conservation": conservation, "dimensional_reduction": reduction,
        "grid_results": grid_results, "grid_comparison_medium_fine": grid_comparison,
        "time_results": time_results, "time_comparison_half_quarter": time_comparison,
    }
    (output / "validation_3d.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    rows = []
    for family, comparison in (("grid", grid_comparison), ("time", time_comparison)):
        for metric, item in comparison.items():
            rows.append({"family": family, "metric": metric, **item})
    with (output / "convergence.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("family", "metric", "first", "second",
            "absolute_error", "relative_error", "excited", "passed"))
        writer.writeheader(); writer.writerows(rows)
    lines = ["# 3D Gate D validation", "", "Evidence: {}".format(EVIDENCE_LABEL), "",
             "Gate D passed: **{}**".format(gate_passed), "",
             "- Conservation: {}".format(conservation_passed),
             "- Dimensional reduction: {}".format(reduction_passed),
             "- Grid/time convergence: {}".format(convergence_passed)]
    (output / "VALIDATION_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report

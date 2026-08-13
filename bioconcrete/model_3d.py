"""Coupled reference 3D reactive-transport solver."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
from pathlib import Path
import time as clock
from typing import Any, Dict, Optional, Tuple

import numpy as np

from .boundary_3d import BoundaryCondition3D
from .chemistry import GeochemLookup, carbonate_fractions
from .config import ModelConfig
from .deposition_3d import DepositionResult3D, update_aperture_3d
from .grid_3d import StructuredGrid3D, capsule_profile_3d, rectangular_grid_3d
from .io_3d import estimate_peak_memory_gb, read_checkpoint, write_checkpoint
from .model import SECONDS_PER_DAY, _initial_state, _is_wet, _state_ph
from .reaction_kernel import reaction_step_cells
from .state import S, STATE_NAMES
from .transport_3d import transport_step_3d


EVIDENCE_LABEL = "uncalibrated 3D model output; not experimental data"
TRANSPORTED = {
    "lactate_mol_m3": "lactate",
    "oxygen_mol_m3": "oxygen",
    "calcium_mol_m3": "calcium",
    "inorganic_carbon_mol_m3": "carbon",
}


@dataclass
class SimulationResult3D:
    times_d: np.ndarray
    state: np.ndarray
    geometry: StructuredGrid3D
    deposition: DepositionResult3D
    aperture_history_m: np.ndarray
    closure_history: np.ndarray
    summary: Dict[str, float]
    diagnostics: Dict[str, object]
    performance: Dict[str, object]
    config: ModelConfig


def inventory_3d(state: np.ndarray, grid: StructuredGrid3D) -> Dict[str, float]:
    volume = grid.cell_volume_m3.ravel()
    cap = state[:, S["capsule_calcium_lactate_mol_m3"]]
    return {
        "carbon_mol": float(np.sum((
            6.0 * cap + 3.0 * state[:, S["lactate_mol_m3"]]
            + state[:, S["inorganic_carbon_mol_m3"]]
            + state[:, S["calcite_mol_m3"]]
            + state[:, S["biomass_carbon_mol_m3"]]
        ) * volume)),
        "calcium_mol": float(np.sum((
            cap + state[:, S["calcium_mol_m3"]]
            + state[:, S["portlandite_mol_m3"]]
            + state[:, S["calcite_mol_m3"]]
        ) * volume)),
        "calcite_mol": float(np.sum(state[:, S["calcite_mol_m3"]] * volume)),
        "csh_volume_m3": float(np.sum(state[:, S["csh_volume_fraction"]] * volume)),
    }


def _new_ledger(state: np.ndarray, grid: StructuredGrid3D) -> Dict[str, Any]:
    inventory = inventory_3d(state, grid)
    return {
        "initial_inventory": inventory,
        "boundary_mol": {species: {face: 0.0 for face in
            ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")}
            for species in TRANSPORTED},
        "reaction_inventory_change": {"carbon_mol": 0.0, "calcium_mol": 0.0},
        "transport_inventory_change": {"carbon_mol": 0.0, "calcium_mol": 0.0},
    }


def _wet_factor(time_s: float, config: ModelConfig) -> float:
    return 1.0 if _is_wet(time_s, config) else config.transport.dry_diffusivity_factor


def _next_environment_event(time_s: float, config: ModelConfig) -> float:
    if config.environment.exposure != "intermittent":
        return float("inf")
    day_start = np.floor(time_s / SECONDS_PER_DAY) * SECONDS_PER_DAY
    wet_end = day_start + config.environment.wet_hours_per_day * 3600.0
    day_end = day_start + SECONDS_PER_DAY
    if time_s < wet_end - 1e-9:
        return wet_end
    return day_end


def _obstruction(state: np.ndarray, grid: StructuredGrid3D, config: ModelConfig) -> Tuple[np.ndarray, DepositionResult3D]:
    deposition = update_aperture_3d(state, grid, config)
    b0 = config.transport.crack_width_mm * 1e-3
    column = np.maximum(deposition.aperture_m / max(b0, 1e-30), 0.0)
    factor = np.broadcast_to(column[:, None, :], grid.shape).copy()
    factor **= config.transport.tortuosity_exponent
    factor[np.broadcast_to(deposition.sealed_mask[:, None, :], grid.shape)] = 0.0
    return factor, deposition


def _boundaries(name: str, config: ModelConfig) -> Dict[str, BoundaryCondition3D]:
    if config.simulation.closed_system:
        return {}
    if name == "oxygen_mol_m3" and config.boundary_3d.oxygen_supply_mode == "boundary_robin":
        return {"x_min": BoundaryCondition3D(
            "robin", config.environment.oxygen_boundary_mol_m3,
            config.boundary_3d.oxygen_mass_transfer_m_s,
        )}
    if name == "inorganic_carbon_mol_m3":
        return {"x_min": BoundaryCondition3D(
            "robin", config.environment.inorganic_carbon_boundary_mol_m3,
            config.boundary_3d.carbon_mass_transfer_m_s,
        )}
    return {}


def _reconcile_carbonate(state: np.ndarray, config: ModelConfig) -> None:
    if config.simulation.carbonate_mode != "equilibrium":
        state[:, S["hydrated_carbon_mol_m3"]] = np.minimum(
            state[:, S["hydrated_carbon_mol_m3"]],
            state[:, S["inorganic_carbon_mol_m3"]],
        )
        return
    ph = _state_ph(state, config)
    _, alpha_hco3, alpha_co3 = carbonate_fractions(ph, config.environment.temperature_c)
    state[:, S["hydrated_carbon_mol_m3"]] = (
        alpha_hco3 + alpha_co3
    ) * state[:, S["inorganic_carbon_mol_m3"]]


def _species_transport(
    state: np.ndarray, grid: StructuredGrid3D, time_s: float, dt_s: float,
    config: ModelConfig, ledger: Dict[str, Any], performance: Dict[str, Any],
) -> np.ndarray:
    base_diffusivity = {
        "lactate_mol_m3": config.transport.diffusivity_lactate_m2_s,
        "oxygen_mol_m3": config.transport.diffusivity_oxygen_m2_s,
        "calcium_mol_m3": config.transport.diffusivity_calcium_m2_s,
        "inorganic_carbon_mol_m3": config.transport.diffusivity_carbon_m2_s,
    }
    before = inventory_3d(state, grid)
    obstruction, _ = _obstruction(state, grid, config)
    wet_factor = _wet_factor(time_s + 0.5 * dt_s, config)
    for name, diffusivity in base_diffusivity.items():
        boundary = _boundaries(name, config)
        field = state[:, S[name]].reshape(grid.shape)
        field, diagnostics = transport_step_3d(
            field, grid, dt_s, diffusivity * wet_factor * obstruction,
            boundary, config.transport.advective_velocity_m_s,
            config.solver_3d.linear_solver, config.solver_3d.relative_tolerance,
            config.solver_3d.absolute_tolerance,
            config.solver_3d.maximum_linear_iterations,
        )
        state[:, S[name]] = field.ravel()
        for face, rate in diagnostics.boundary_rates_after.items():
            ledger["boundary_mol"][name][face] += rate * dt_s
        performance["linear_iterations"] += diagnostics.linear.iterations
        performance["linear_solves"] += 1
    _reconcile_carbonate(state, config)
    after = inventory_3d(state, grid)
    ledger["transport_inventory_change"]["carbon_mol"] += after["carbon_mol"] - before["carbon_mol"]
    ledger["transport_inventory_change"]["calcium_mol"] += after["calcium_mol"] - before["calcium_mol"]
    return state


def _reaction(
    state: np.ndarray, grid: StructuredGrid3D, time_s: float, dt_s: float,
    config: ModelConfig, geochem: Optional[GeochemLookup], ledger: Dict[str, Any],
) -> np.ndarray:
    before = inventory_3d(state, grid)
    state = reaction_step_cells(
        state, time_s, dt_s, config, geochem,
        batch_size=config.solver_3d.reaction_batch_cells,
        workers=config.solver_3d.reaction_workers,
        parallel_backend=config.solver_3d.reaction_parallel_backend,
    )
    after = inventory_3d(state, grid)
    ledger["reaction_inventory_change"]["carbon_mol"] += after["carbon_mol"] - before["carbon_mol"]
    ledger["reaction_inventory_change"]["calcium_mol"] += after["calcium_mol"] - before["calcium_mol"]
    return state


def _closure(ledger: Dict[str, Any], final: Dict[str, float]) -> Dict[str, Any]:
    carbon_boundary = (
        3.0 * sum(ledger["boundary_mol"]["lactate_mol_m3"].values())
        + sum(ledger["boundary_mol"]["inorganic_carbon_mol_m3"].values())
    )
    calcium_boundary = sum(ledger["boundary_mol"]["calcium_mol_m3"].values())
    result = {}
    for key, boundary in (("carbon_mol", carbon_boundary), ("calcium_mol", calcium_boundary)):
        initial = ledger["initial_inventory"][key]
        change = final[key] - initial
        reaction = ledger["reaction_inventory_change"][key]
        residual = change - boundary - reaction
        result[key] = {
            "initial": initial, "final": final[key], "state_change": change,
            "boundary_integral": boundary, "reaction_integral": reaction,
            "residual": residual,
            "relative_error": abs(residual) / max(abs(initial), abs(change), 1e-30),
        }
    return result


def simulate_3d(
    config: Optional[ModelConfig] = None,
    geochem: Optional[GeochemLookup] = None,
    geometry: Optional[StructuredGrid3D] = None,
    checkpoint: Optional[Path] = None,
    resume: bool = False,
    failure_manifest: Optional[Path] = None,
    stop_after_s: Optional[float] = None,
) -> SimulationResult3D:
    started = clock.perf_counter()
    config = deepcopy(config or ModelConfig())
    config.validate()
    grid = geometry or rectangular_grid_3d(config)
    saved_count = int(np.ceil(config.simulation.days / config.output_3d.save_every_days)) + 1
    memory_gb = estimate_peak_memory_gb(
        grid, saved_count, len(STATE_NAMES), config.output_3d.save_full_state
    )
    if memory_gb > config.solver_3d.memory_limit_gb:
        raise MemoryError("estimated {:.3f} GB exceeds solver_3d.memory_limit_gb".format(memory_gb))

    profile, _ = capsule_profile_3d(config, grid)
    state = _initial_state(config, profile.ravel(), "3d")
    current_time = 0.0
    step = output_index = retries = 0
    ledger = _new_ledger(state, grid)
    if resume:
        if checkpoint is None:
            raise ValueError("resume requires checkpoint")
        current_time, state, metadata = read_checkpoint(checkpoint, config, grid, full=True)
        ledger = metadata["ledger"]
        step = int(metadata["step"])
        output_index = int(metadata["output_index"])
        retries = int(metadata["retry_count"])

    reaction_config = deepcopy(config)
    if config.boundary_3d.oxygen_supply_mode == "boundary_robin":
        reaction_config.environment.oxygen_transfer_s = 0.0
    configured_total_s = config.simulation.days * SECONDS_PER_DAY
    total_s = min(configured_total_s, stop_after_s) if stop_after_s is not None else configured_total_s
    if total_s < current_time:
        raise ValueError("stop_after_s precedes checkpoint time")
    output_interval = max(config.output_3d.save_every_days * SECONDS_PER_DAY, 1.0)
    targets = np.unique(np.append(np.arange(current_time, total_s, output_interval), total_s))
    snapshots = [state.copy()]
    aperture = update_aperture_3d(state, grid, config)
    aperture_history = [aperture.aperture_m.copy()]
    closure_history = [aperture.closure.copy()]
    performance = {"linear_iterations": 0, "linear_solves": 0, "time_step_retries": retries,
                   "checkpoint_count": 0, "reaction_time_s": 0.0,
                   "transport_time_s": 0.0}
    maximum_step = config.simulation.reaction_step_h * 3600.0
    checkpoint_stride = max(config.solver_3d.checkpoint_interval_steps, 1)
    max_retries = 4

    for target in targets[1:]:
        while current_time < target - 1e-9:
            event = min(target, _next_environment_event(current_time, config))
            attempted = min(maximum_step, event - current_time)
            retry = 0
            while True:
                candidate = state.copy()
                candidate_ledger = deepcopy(ledger)
                try:
                    if config.solver_3d.splitting_scheme == "strang":
                        stamp = clock.perf_counter()
                        candidate = _species_transport(candidate, grid, current_time, attempted / 2,
                                                       config, candidate_ledger, performance)
                        performance["transport_time_s"] += clock.perf_counter() - stamp
                    stamp = clock.perf_counter()
                    candidate = _reaction(candidate, grid, current_time, attempted,
                                          reaction_config, geochem, candidate_ledger)
                    performance["reaction_time_s"] += clock.perf_counter() - stamp
                    stamp = clock.perf_counter()
                    candidate = _species_transport(
                        candidate, grid, current_time + attempted / 2,
                        attempted / 2 if config.solver_3d.splitting_scheme == "strang" else attempted,
                        config, candidate_ledger, performance,
                    )
                    performance["transport_time_s"] += clock.perf_counter() - stamp
                    if not np.isfinite(candidate).all() or np.min(candidate) < 0:
                        raise RuntimeError("non-finite or negative coupled state")
                    state, ledger = candidate, candidate_ledger
                    break
                except Exception as error:
                    retry += 1
                    retries += 1
                    performance["time_step_retries"] = retries
                    if retry > max_retries:
                        path = failure_manifest or Path("model_runs/v0.6.0/3d/failure_manifest.json")
                        path.parent.mkdir(parents=True, exist_ok=True)
                        path.write_text(json.dumps({
                            "status": "failed", "time_s": current_time,
                            "attempted_dt_s": attempted, "retries": retry,
                            "error": repr(error), "evidence_label": EVIDENCE_LABEL,
                        }, indent=2), encoding="utf-8")
                        raise RuntimeError("3D time step failed; see {}".format(path)) from error
                    attempted *= 0.5
            current_time += attempted
            step += 1
            if checkpoint is not None and step % checkpoint_stride == 0:
                write_checkpoint(checkpoint, current_time, state, config, grid, ledger,
                                 step, output_index, retries)
                performance["checkpoint_count"] += 1
        output_index += 1
        snapshots.append(state.copy())
        deposition = update_aperture_3d(state, grid, config)
        aperture_history.append(deposition.aperture_m.copy())
        closure_history.append(deposition.closure.copy())

    if checkpoint is not None:
        write_checkpoint(checkpoint, current_time, state, config, grid, ledger,
                         step, output_index, retries)
        performance["checkpoint_count"] += 1
    deposition = update_aperture_3d(state, grid, config)
    final_inventory = inventory_3d(state, grid)
    balance = _closure(ledger, final_inventory)
    oxygen = state[:, S["oxygen_mol_m3"]].reshape(grid.shape)
    summary = {
        "calcite_mol": final_inventory["calcite_mol"],
        "csh_volume_m3": final_inventory["csh_volume_m3"],
        "area_weighted_closure": deposition.area_weighted_closure,
        "maximum_local_closure": float(np.max(deposition.closure)),
        "open_volume_m3": float(np.sum(deposition.aperture_m *
            (grid.face_area_y_m2[:, 0, :]))),
        "open_volume_closure": deposition.open_volume_closure,
        "sealed_area_fraction": deposition.sealed_area_fraction,
        "minimum_aperture_m": float(np.min(deposition.aperture_m)),
        "oxygen_penetration_depth_m": float(grid.x_m[np.max(np.where(
            oxygen.mean(axis=(0, 1)) > 0.01 * max(config.environment.oxygen_boundary_mol_m3, 1e-30)
        )[0])] if np.any(oxygen.mean(axis=(0, 1)) > 0.01 * max(config.environment.oxygen_boundary_mol_m3, 1e-30)) else 0.0),
    }
    diagnostics = {
        "axis_order": ["z", "y", "x"], "flatten_order": "C",
        "geometry_hash": grid.geometry_hash, "evidence_label": EVIDENCE_LABEL,
        "nonnegative": bool(np.min(state) >= 0), "finite": bool(np.isfinite(state).all()),
        "z_oxygen_range_mol_m3": float(np.ptp(oxygen.mean(axis=(1, 2)))),
        "estimated_peak_memory_gb": memory_gb, "steps": step,
        "conservation": balance, "ledger": ledger,
    }
    performance["wall_time_s"] = clock.perf_counter() - started
    return SimulationResult3D(
        targets / SECONDS_PER_DAY,
        np.stack(snapshots).reshape((len(targets),) + grid.shape + (len(STATE_NAMES),)),
        grid, deposition, np.stack(aperture_history), np.stack(closure_history),
        summary, diagnostics, performance, config,
    )

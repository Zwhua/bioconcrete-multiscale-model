"""0D, 1D, and 2D coupled reaction-transport solvers."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import spsolve

from .chemistry import GeochemLookup, carbonate_fractions, ph_from_alkalinity
from .config import ModelConfig


SECONDS_PER_DAY = 86400.0
STATE_NAMES = (
    "capsule_calcium_lactate_mol_m3",
    "spore_density_rel",
    "active_density_rel",
    "lactate_mol_m3",
    "oxygen_mol_m3",
    "calcium_mol_m3",
    "inorganic_carbon_mol_m3",
    "hydrated_carbon_mol_m3",
    "portlandite_mol_m3",
    "calcite_mol_m3",
    "csh_volume_fraction",
    "biomass_carbon_mol_m3",
    "ammonium_mol_m3",
    "total_alkalinity_mol_m3",
    "environment_signal",
    "activation_state",
    "activation_memory_h",
    "tracked_oxygen_mol_m3",
    "tracked_ph",
    "cumulative_activity_h",
    "premature_consumption_mol_m3",
    "activation_delay_h",
)
S = {name: index for index, name in enumerate(STATE_NAMES)}
DISSOLVED = {
    S["lactate_mol_m3"]: "lactate",
    S["oxygen_mol_m3"]: "oxygen",
    S["calcium_mol_m3"]: "calcium",
    S["inorganic_carbon_mol_m3"]: "carbon",
    S["hydrated_carbon_mol_m3"]: "carbon",
}


@dataclass
class SimulationResult:
    level: str
    frame: pd.DataFrame
    summary: Dict[str, float]
    diagnostics: Dict[str, object]
    config: ModelConfig

    def save(self, output_dir: Path) -> None:
        output_dir.mkdir(parents=True, exist_ok=True)
        self.frame.to_csv(output_dir / "state.csv", index=False)
        (output_dir / "summary.json").write_text(json.dumps(self.summary, indent=2), encoding="utf-8")
        (output_dir / "diagnostics.json").write_text(json.dumps(self.diagnostics, indent=2), encoding="utf-8")
        self.config.save(output_dir / "config.json")


def _sigmoid(value: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(value, -60.0, 60.0)))


def _is_wet(time_s: float, config: ModelConfig) -> bool:
    exposure = config.environment.exposure
    if exposure == "continuous":
        return True
    if exposure == "dry":
        return False
    hour = (time_s / 3600.0) % 24.0
    return hour < config.environment.wet_hours_per_day


def _environment(time_s: float, config: ModelConfig) -> Tuple[float, float]:
    wet = _is_wet(time_s, config)
    aw = config.environment.water_activity_wet if wet else config.environment.water_activity_dry
    return aw, 1.0 if wet else config.transport.dry_diffusivity_factor


def _state_ph(state: np.ndarray, config: ModelConfig) -> np.ndarray:
    """Return fixed or charge-balanced pH for each reaction cell."""

    if config.simulation.ph_mode == "fixed":
        return np.full(state.shape[0], config.environment.ph, dtype=float)
    return ph_from_alkalinity(
        state[:, S["inorganic_carbon_mol_m3"]],
        state[:, S["total_alkalinity_mol_m3"]],
        config.environment.temperature_c,
        config.environment.ph_minimum,
        config.environment.ph_maximum,
        strict=False,
    )


def _environment_suitability(
    aw: float, oxygen: np.ndarray, ph: np.ndarray, config: ModelConfig
) -> np.ndarray:
    k = config.kinetics
    env = config.environment
    water_gate = _sigmoid(80.0 * (aw - k.aw_threshold))
    oxygen_gate = _sigmoid(35.0 * (oxygen - k.oxygen_threshold_mol_m3))
    ph_gate = np.exp(-0.5 * ((ph - k.ph_optimum) / max(k.ph_width, 1e-6)) ** 2)
    temperature_gate = np.exp(
        -0.5 * ((env.temperature_c - k.temperature_optimum_c) / max(k.temperature_width_c, 1e-6)) ** 2
    )
    return water_gate * oxygen_gate * ph_gate * temperature_gate


def _effective_dosage(config: ModelConfig, level: str) -> float:
    """Return concentration scaling for the configured material-dose basis."""

    dosage = config.kinetics.agent_dosage_multiplier
    if config.kinetics.dosage_basis == "fixed_total_inventory":
        dosage *= (
            config.kinetics.reference_inventory_volume_m3
            / _total_crack_volume_m3(config, level)
        )
    return dosage


def _initial_state(
    config: ModelConfig, capsule_profile: np.ndarray, level: str = "0d"
) -> np.ndarray:
    n_cells = capsule_profile.size
    state = np.zeros((n_cells, len(STATE_NAMES)), dtype=float)
    dosage = _effective_dosage(config, level)
    state[:, S["capsule_calcium_lactate_mol_m3"]] = (
        config.kinetics.capsule_calcium_lactate_mol_m3 * dosage * capsule_profile
    )
    state[:, S["spore_density_rel"]] = config.kinetics.spore_density_rel * dosage * capsule_profile
    state[:, S["active_density_rel"]] = config.kinetics.active_density_rel * dosage * capsule_profile
    state[:, S["oxygen_mol_m3"]] = config.environment.oxygen_initial_mol_m3
    state[:, S["portlandite_mol_m3"]] = config.chemistry.portlandite_mol_m3
    state[:, S["total_alkalinity_mol_m3"]] = config.environment.initial_alkalinity_mol_m3
    state[:, S["tracked_oxygen_mol_m3"]] = 0.0
    state[:, S["tracked_ph"]] = min(config.environment.ph_maximum, config.environment.ph + 0.5)
    return state


def _reaction_rhs(
    _time: float,
    flat_state: np.ndarray,
    n_cells: int,
    config: ModelConfig,
    aw: float,
    wet: bool,
    geochem: Optional[GeochemLookup],
) -> np.ndarray:
    y = np.nan_to_num(
        np.maximum(flat_state.reshape(n_cells, len(STATE_NAMES)), 0.0),
        nan=0.0,
        posinf=1e8,
        neginf=0.0,
    )
    dy = np.zeros_like(y)
    k = config.kinetics
    chem = config.chemistry
    env = config.environment

    capsule = y[:, S["capsule_calcium_lactate_mol_m3"]]
    spores = y[:, S["spore_density_rel"]]
    active = y[:, S["active_density_rel"]]
    lactate = y[:, S["lactate_mol_m3"]]
    oxygen = y[:, S["oxygen_mol_m3"]]
    calcium = y[:, S["calcium_mol_m3"]]
    carbon = y[:, S["inorganic_carbon_mol_m3"]]
    ph = _state_ph(y, config)
    _, alpha_hco3, alpha_co3 = carbonate_fractions(ph, env.temperature_c)
    hydrated_target_fraction = alpha_hco3 + alpha_co3
    if config.simulation.carbonate_mode == "equilibrium":
        hydrated = hydrated_target_fraction * carbon
    else:
        hydrated = np.minimum(y[:, S["hydrated_carbon_mol_m3"]], carbon)
    portlandite = y[:, S["portlandite_mol_m3"]]

    signal = y[:, S["environment_signal"]]
    activation = y[:, S["activation_state"]]
    memory_h = y[:, S["activation_memory_h"]]
    tracked_oxygen = y[:, S["tracked_oxygen_mol_m3"]]
    tracked_ph = y[:, S["tracked_ph"]]
    tracking_tau_s = max(k.signal_relaxation_h * 3600.0, 1.0)
    oxygen_rise_h = (oxygen - tracked_oxygen) / tracking_tau_s * 3600.0
    ph_drop_h = (tracked_ph - ph) / tracking_tau_s * 3600.0
    oxygen_change_gate = _sigmoid(
        80.0 * (oxygen_rise_h - k.oxygen_rise_threshold_mol_m3_h)
    )
    ph_change_gate = _sigmoid(8.0 * (ph_drop_h - k.ph_drop_threshold_h))
    suitability = _environment_suitability(aw, oxygen, ph, config)
    if k.gate_logic == "AND":
        signal_target = suitability * oxygen_change_gate * ph_change_gate
    elif k.gate_logic == "OR":
        change_gate = 1.0 - (1.0 - oxygen_change_gate) * (1.0 - ph_change_gate)
        signal_target = suitability * change_gate
    else:
        signal_target = suitability
    duration_scale_h = max(k.activation_duration_h, 1.0e-6)
    memory_target = _sigmoid(
        8.0 * (memory_h - k.activation_duration_h) / duration_scale_h
    )
    response_tau_s = max(k.response_delay_h * 3600.0, 1.0)
    effective_gate = k.basal_leak_fraction + (1.0 - k.basal_leak_fraction) * np.clip(activation, 0.0, 1.0)
    release_gate = _sigmoid(np.full(n_cells, 80.0 * (aw - k.aw_threshold)))
    release = k.capsule_release_s * release_gate * capsule
    germination = k.germination_s * effective_gate * spores

    monod_l = lactate / (k.effective_km_mol_m3 + lactate + 1e-30)
    monod_o = oxygen / (k.k_oxygen_mol_m3 + oxygen + 1e-30)
    effective_activity = (
        k.effective_kcat_s
        * k.active_unit_concentration
        * k.activity_multiplier
        * active
    )
    uptake = effective_activity * monod_l * monod_o * effective_gate
    carbon_fraction = k.biomass_carbon_fraction
    aerobic_fraction = max(1.0 - carbon_fraction, 1e-9)
    uptake = np.minimum(uptake, lactate / 300.0)
    uptake = np.minimum(uptake, oxygen / (3.0 * aerobic_fraction * 300.0 + 1e-30))

    if config.simulation.carbonate_mode == "equilibrium":
        hydration = np.zeros(n_cells)
        dehydration = np.zeros(n_cells)
    else:
        hydration = chem.ca_hydration_s * np.maximum(hydrated_target_fraction * carbon - hydrated, 0.0)
        dehydration = chem.ca_dehydration_s * np.maximum(hydrated - hydrated_target_fraction * carbon, 0.0)

    carbonate_available = hydrated * alpha_co3 / np.maximum(hydrated_target_fraction, 1e-30)
    lookup = geochem or GeochemLookup()
    omega = np.clip(
        np.nan_to_num(lookup.saturation(calcium, carbon, ph, chem, env.temperature_c), nan=0.0),
        0.0,
        1e12,
    )
    precipitation = (
        chem.calcite_rate_mol_m3_s
        * chem.calcite_surface_area_rel
        * np.power(np.maximum(omega - 1.0, 0.0), chem.calcite_reaction_order)
    )
    precipitation = np.minimum.reduce(
        [precipitation, calcium / 300.0, carbon / 300.0, carbonate_available / 300.0]
    )
    ch_dissolution = (
        chem.portlandite_dissolution_s
        * portlandite
        * carbon
        / (chem.portlandite_carbon_half_mol_m3 + carbon + 1e-30)
    )
    ch_dissolution = np.minimum(ch_dissolution, portlandite / 300.0)

    alkaline_excess = np.maximum(ph - 12.0, 0.0)
    decay_rate = k.decay_s + k.alkaline_decay_s * alkaline_excess
    encapsulation_decay = k.encapsulation_decay_m3_mol * precipitation * active
    growth = k.maximum_growth_s * monod_l * monod_o * effective_gate * active + k.biomass_yield_rel_m3_mol * uptake

    dy[:, S["capsule_calcium_lactate_mol_m3"]] = -release
    dy[:, S["spore_density_rel"]] = -germination - decay_rate * spores
    dy[:, S["active_density_rel"]] = germination + growth - decay_rate * active - encapsulation_decay
    dy[:, S["lactate_mol_m3"]] = 2.0 * release - uptake
    dy[:, S["oxygen_mol_m3"]] = -3.0 * aerobic_fraction * uptake
    if wet and not config.simulation.closed_system:
        dy[:, S["oxygen_mol_m3"]] += env.oxygen_transfer_s * np.maximum(
            env.oxygen_boundary_mol_m3 - oxygen, 0.0
        )
    dy[:, S["calcium_mol_m3"]] = release + ch_dissolution - precipitation
    carbon_change = 3.0 * aerobic_fraction * uptake - precipitation
    dy[:, S["inorganic_carbon_mol_m3"]] = carbon_change
    if config.simulation.carbonate_mode == "equilibrium":
        dy[:, S["hydrated_carbon_mol_m3"]] = hydrated_target_fraction * carbon_change
    else:
        dy[:, S["hydrated_carbon_mol_m3"]] = hydration - dehydration - precipitation
    dy[:, S["portlandite_mol_m3"]] = -ch_dissolution
    dy[:, S["calcite_mol_m3"]] = precipitation
    dy[:, S["csh_volume_fraction"]] = (
        k.csh_release_s
        * release_gate
        * capsule
        / max(k.capsule_calcium_lactate_mol_m3, 1e-30)
        * k.capsule_csh_volume_fraction
    )
    dy[:, S["biomass_carbon_mol_m3"]] = 3.0 * carbon_fraction * uptake
    # The selected pathway is ammonia-free. This state is retained as an invariant diagnostic.
    dy[:, S["ammonium_mol_m3"]] = 0.0
    if config.simulation.ph_mode == "dynamic":
        # Portlandite dissolution contributes two charge equivalents; calcite
        # precipitation consumes two. This is the model's conserved alkalinity state.
        dy[:, S["total_alkalinity_mol_m3"]] = 2.0 * ch_dissolution - 2.0 * precipitation
    dy[:, S["environment_signal"]] = (signal_target - signal) / tracking_tau_s
    accumulating = np.clip(signal, 0.0, 1.0) / 3600.0
    relaxing = np.maximum(memory_h, 0.0) / max(k.activation_duration_h * 3600.0, 1.0)
    dy[:, S["activation_memory_h"]] = np.where(signal >= 0.5, accumulating, -relaxing)
    dy[:, S["activation_state"]] = (memory_target - activation) / response_tau_s
    dy[:, S["tracked_oxygen_mol_m3"]] = (oxygen - tracked_oxygen) / tracking_tau_s
    dy[:, S["tracked_ph"]] = (ph - tracked_ph) / tracking_tau_s
    dy[:, S["cumulative_activity_h"]] = effective_gate / 3600.0
    dy[:, S["premature_consumption_mol_m3"]] = uptake * (1.0 - signal_target)
    dy[:, S["activation_delay_h"]] = (
        ((signal_target >= 0.5) & (activation < 0.5)).astype(float) / 3600.0
    )
    return dy.ravel()


def _reaction_step(
    state: np.ndarray,
    time_s: float,
    dt_s: float,
    config: ModelConfig,
    geochem: Optional[GeochemLookup],
) -> np.ndarray:
    n_cells = state.shape[0]
    aw, _ = _environment(time_s + 0.5 * dt_s, config)
    wet = _is_wet(time_s + 0.5 * dt_s, config)
    options = {}
    if n_cells > 4:
        block = sparse.csr_matrix(np.ones((len(STATE_NAMES), len(STATE_NAMES)), dtype=bool))
        options["jac_sparsity"] = sparse.kron(sparse.eye(n_cells, format="csr"), block, format="csr")
    try:
        solution = solve_ivp(
            _reaction_rhs,
            (0.0, dt_s),
            state.ravel(),
            method="BDF",
            rtol=2.0e-5,
            atol=1.0e-9,
            args=(n_cells, config, aw, wet, geochem),
            **options
        )
        if not solution.success or not np.isfinite(solution.y[:, -1]).all():
            raise RuntimeError(solution.message)
        updated = np.maximum(solution.y[:, -1].reshape(state.shape), 0.0)
    except (RuntimeError, ValueError):
        # Local reactions are independent during the split step. Cell-wise BDF is
        # slower but avoids a single ill-conditioned block stopping a 2D run.
        updated = np.empty_like(state)
        for cell in range(n_cells):
            local = solve_ivp(
                _reaction_rhs,
                (0.0, dt_s),
                state[cell],
                method="BDF",
                rtol=2.0e-5,
                atol=1.0e-9,
                args=(1, config, aw, wet, geochem),
            )
            if not local.success:
                raise RuntimeError("Reaction solve failed in cell {}: {}".format(cell, local.message))
            updated[cell] = np.maximum(local.y[:, -1], 0.0)
    if config.simulation.carbonate_mode == "equilibrium":
        ph = _state_ph(updated, config)
        _, alpha_hco3, alpha_co3 = carbonate_fractions(
            ph, config.environment.temperature_c
        )
        updated[:, S["hydrated_carbon_mol_m3"]] = (alpha_hco3 + alpha_co3) * updated[
            :, S["inorganic_carbon_mol_m3"]
        ]
    else:
        updated[:, S["hydrated_carbon_mol_m3"]] = np.minimum(
            updated[:, S["hydrated_carbon_mol_m3"]], updated[:, S["inorganic_carbon_mol_m3"]]
        )
    updated[:, S["ammonium_mol_m3"]] = 0.0
    updated[:, S["environment_signal"]] = np.clip(updated[:, S["environment_signal"]], 0.0, 1.0)
    updated[:, S["activation_state"]] = np.clip(updated[:, S["activation_state"]], 0.0, 1.0)
    if config.simulation.ph_mode == "dynamic":
        # Accepted steps must satisfy the configured charge-balance interval.
        ph_from_alkalinity(
            updated[:, S["inorganic_carbon_mol_m3"]],
            updated[:, S["total_alkalinity_mol_m3"]],
            config.environment.temperature_c,
            config.environment.ph_minimum,
            config.environment.ph_maximum,
            strict=True,
        )
    return updated


def _uniform_cell_geometry(config: ModelConfig, level: str, n_cells: int) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-cell crack volume and total opposing-wall area."""

    total_volume = _total_crack_volume_m3(config, level)
    aperture_m = config.transport.crack_width_mm * 1.0e-3
    cell_volume = np.full(n_cells, total_volume / max(n_cells, 1))
    total_wall_area = 2.0 * cell_volume / aperture_m
    return cell_volume, total_wall_area


def repair_metrics(
    state: np.ndarray,
    config: ModelConfig,
    cell_volume_m3: Optional[np.ndarray] = None,
    total_wall_area_m2: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """Map solid volume to two-wall deposition and aperture closure.

    ``total_wall_area_m2`` is the sum of both opposing crack-wall areas.
    ``wall_deposition_thickness_mm`` is the thickness on one wall, so closure
    uses twice that thickness exactly once.
    """

    chem = config.chemistry
    trans = config.transport
    calcite_volume = state[:, S["calcite_mol_m3"]] * chem.calcite_molar_mass_kg_mol / chem.calcite_density_kg_m3
    csh_volume = state[:, S["csh_volume_fraction"]]
    solid_fraction = np.clip(calcite_volume + csh_volume, 0.0, 1.0)
    n_cells = state.shape[0]
    if cell_volume_m3 is None:
        cell_volume_m3 = np.full(n_cells, _total_crack_volume_m3(config, "0d"))
    else:
        cell_volume_m3 = np.broadcast_to(np.asarray(cell_volume_m3, dtype=float), (n_cells,))
    if total_wall_area_m2 is None:
        total_wall_area_m2 = np.full(
            n_cells, 2.0 * trans.crack_length_mm * trans.crack_depth_mm * 1.0e-6
        )
    else:
        total_wall_area_m2 = np.broadcast_to(np.asarray(total_wall_area_m2, dtype=float), (n_cells,))
    calcite_solid_volume_m3 = calcite_volume * cell_volume_m3
    csh_solid_volume_m3 = csh_volume * cell_volume_m3
    total_solid_volume_m3 = solid_fraction * cell_volume_m3
    calcite_wall_volume_m3 = calcite_solid_volume_m3 * chem.wall_deposition_fraction
    csh_wall_volume_m3 = csh_solid_volume_m3 * chem.wall_deposition_fraction
    wall_solid_volume_m3 = total_solid_volume_m3 * chem.wall_deposition_fraction
    nonwall_solid_volume_m3 = total_solid_volume_m3 - wall_solid_volume_m3
    one_wall_thickness_m = wall_solid_volume_m3 / np.maximum(total_wall_area_m2, 1e-30)
    wall_deposition_thickness_mm = one_wall_thickness_m * 1.0e3
    calcite_closure_contribution = np.clip(
        2.0 * calcite_wall_volume_m3
        / np.maximum(total_wall_area_m2, 1e-30)
        * 1.0e3 / trans.crack_width_mm,
        0.0, 1.0,
    )
    csh_closure_contribution = np.clip(
        2.0 * csh_wall_volume_m3
        / np.maximum(total_wall_area_m2, 1e-30)
        * 1.0e3 / trans.crack_width_mm,
        0.0, 1.0,
    )
    crack_closure_ratio = np.clip(
        2.0 * wall_deposition_thickness_mm / trans.crack_width_mm, 0.0, 1.0
    )
    aperture = trans.crack_width_mm * (1.0 - crack_closure_ratio)
    porosity = np.maximum(
        trans.porosity_initial - solid_fraction * trans.porosity_initial,
        trans.porosity_minimum,
    )
    phi0 = trans.porosity_initial
    permeability = (porosity / phi0) ** 3 * ((1.0 - phi0) / np.maximum(1.0 - porosity, 1e-12)) ** 2
    transmissivity = np.maximum(1.0 - crack_closure_ratio, 0.0) ** 3
    return {
        "solid_fill_fraction": solid_fraction,
        "total_solid_volume_m3": total_solid_volume_m3,
        "calcite_solid_volume_m3": calcite_solid_volume_m3,
        "csh_solid_volume_m3": csh_solid_volume_m3,
        "wall_solid_volume_m3": wall_solid_volume_m3,
        "calcite_wall_volume_m3": calcite_wall_volume_m3,
        "csh_wall_volume_m3": csh_wall_volume_m3,
        "nonwall_solid_volume_m3": nonwall_solid_volume_m3,
        "wall_deposition_thickness_mm": wall_deposition_thickness_mm,
        "calcite_closure_contribution": calcite_closure_contribution,
        "csh_closure_contribution": csh_closure_contribution,
        "crack_closure_ratio": crack_closure_ratio,
        # Deprecated compatibility alias. New reports use crack_closure_ratio.
        "healing_ratio": crack_closure_ratio,
        "aperture_mm": aperture,
        "porosity": porosity,
        "permeability_ratio": permeability,
        "transmissivity_ratio": transmissivity,
        "sorptivity_ratio": np.sqrt(np.maximum(permeability, 0.0)),
    }


def _balance(state: np.ndarray) -> Dict[str, float]:
    cap = state[:, S["capsule_calcium_lactate_mol_m3"]]
    carbon = (
        6.0 * cap
        + 3.0 * state[:, S["lactate_mol_m3"]]
        + state[:, S["inorganic_carbon_mol_m3"]]
        + state[:, S["calcite_mol_m3"]]
        + state[:, S["biomass_carbon_mol_m3"]]
    )
    calcium = cap + state[:, S["calcium_mol_m3"]] + state[:, S["portlandite_mol_m3"]] + state[:, S["calcite_mol_m3"]]
    return {"carbon": float(np.mean(carbon)), "calcium": float(np.mean(calcium))}


def _total_crack_volume_m3(config: ModelConfig, level: str) -> float:
    trans = config.transport
    unresolved_mm = (
        trans.out_of_plane_thickness_mm if level == "2d" else trans.crack_depth_mm
    )
    return (
        trans.crack_length_mm
        * trans.crack_width_mm
        * unresolved_mm
        * 1.0e-9
    )


def _summary(state: np.ndarray, config: ModelConfig, level: str) -> Dict[str, float]:
    cell_volume, wall_area = _uniform_cell_geometry(config, level, state.shape[0])
    metrics = repair_metrics(state, config, cell_volume, wall_area)
    calcite_mean = float(np.mean(state[:, S["calcite_mol_m3"]]))
    return {
        "mean_crack_closure_ratio": float(np.mean(metrics["crack_closure_ratio"])),
        "max_crack_closure_ratio": float(np.max(metrics["crack_closure_ratio"])),
        "mean_healing_ratio": float(np.mean(metrics["crack_closure_ratio"])),
        "mean_calcite_closure_contribution": float(np.mean(metrics["calcite_closure_contribution"])),
        "mean_csh_closure_contribution": float(np.mean(metrics["csh_closure_contribution"])),
        "mean_permeability_ratio": float(np.mean(metrics["permeability_ratio"])),
        "mean_transmissivity_ratio": float(np.mean(metrics["transmissivity_ratio"])),
        "calcite_mol_m3_mean": calcite_mean,
        "calcite_kg_m3_mean": float(
            calcite_mean * config.chemistry.calcite_molar_mass_kg_mol
        ),
        "calcite_mass_mg": calcite_mean
        * _total_crack_volume_m3(config, level)
        * config.chemistry.calcite_molar_mass_kg_mol
        * 1.0e6,
        "ammonium_mol_m3_max": float(np.max(state[:, S["ammonium_mol_m3"]])),
    }


def _state_frame(state: np.ndarray, time_d: float, config: ModelConfig, coordinates: Dict[str, np.ndarray]) -> pd.DataFrame:
    values = {name: state[:, index] for index, name in enumerate(STATE_NAMES)}
    values.update(coordinates)
    n_cells = state.shape[0]
    trans = config.transport
    if "y_mm" in coordinates:
        level = "2d"
    elif "x_mm" in coordinates:
        level = "1d"
    else:
        level = "0d"
    cell_volume, wall_area = _uniform_cell_geometry(config, level, n_cells)
    values.update(repair_metrics(state, config, cell_volume, wall_area))
    values["cell_volume_m3"] = cell_volume
    values["wall_area_m2"] = wall_area
    values["ph"] = _state_ph(state, config)
    effective_gate = config.kinetics.basal_leak_fraction + (
        1.0 - config.kinetics.basal_leak_fraction
    ) * state[:, S["activation_state"]]
    oxygen = state[:, S["oxygen_mol_m3"]]
    ph = values["ph"]
    tracking_tau_s = max(config.kinetics.signal_relaxation_h * 3600.0, 1.0)
    oxygen_rise_h = (
        oxygen - state[:, S["tracked_oxygen_mol_m3"]]
    ) / tracking_tau_s * 3600.0
    ph_drop_h = (
        state[:, S["tracked_ph"]] - ph
    ) / tracking_tau_s * 3600.0
    strict_reference = (
        _environment_suitability(
            config.environment.water_activity_wet if _is_wet(time_d * SECONDS_PER_DAY, config)
            else config.environment.water_activity_dry,
            oxygen, ph, config,
        )
        * _sigmoid(80.0 * (oxygen_rise_h - config.kinetics.oxygen_rise_threshold_mol_m3_h))
        * _sigmoid(8.0 * (ph_drop_h - config.kinetics.ph_drop_threshold_h))
    )
    values["true_activation_index"] = effective_gate * strict_reference
    values["false_activation_index"] = effective_gate * (1.0 - strict_reference)
    values["time_d"] = np.full(state.shape[0], time_d)
    return pd.DataFrame(values)


def _inventory(state: np.ndarray, config: ModelConfig, level: str) -> Dict[str, float]:
    cell_volume, _ = _uniform_cell_geometry(config, level, state.shape[0])
    cap = state[:, S["capsule_calcium_lactate_mol_m3"]]
    carbon = (
        6.0 * cap + 3.0 * state[:, S["lactate_mol_m3"]]
        + state[:, S["inorganic_carbon_mol_m3"]]
        + state[:, S["calcite_mol_m3"]]
        + state[:, S["biomass_carbon_mol_m3"]]
    )
    calcium = (
        cap + state[:, S["calcium_mol_m3"]]
        + state[:, S["portlandite_mol_m3"]]
        + state[:, S["calcite_mol_m3"]]
    )
    return {
        "capsule_calcium_lactate_mol": float(np.sum(cap * cell_volume)),
        "spore_inventory_rel_m3": float(np.sum(state[:, S["spore_density_rel"]] * cell_volume)),
        "active_inventory_rel_m3": float(np.sum(state[:, S["active_density_rel"]] * cell_volume)),
        "remaining_csh_payload_m3": float(np.sum(
            cap / max(config.kinetics.capsule_calcium_lactate_mol_m3, 1e-30)
            * config.kinetics.capsule_csh_volume_fraction * cell_volume
        )),
        "released_csh_m3": float(np.sum(state[:, S["csh_volume_fraction"]] * cell_volume)),
        "carbon_mol": float(np.sum(carbon * cell_volume)),
        "calcium_mol": float(np.sum(calcium * cell_volume)),
        "calcite_mol": float(np.sum(state[:, S["calcite_mol_m3"]] * cell_volume)),
    }


def _diagnostics(
    initial: np.ndarray, final: np.ndarray, config: ModelConfig, level: str
) -> Dict[str, object]:
    before = _balance(initial)
    after = _balance(final)
    errors = {
        key: abs(after[key] - before[key]) / max(abs(before[key]), 1e-30) for key in before
    }
    return {
        "initial_balance": before,
        "final_balance": after,
        "initial_inventory": _inventory(initial, config, level),
        "final_inventory": _inventory(final, config, level),
        "relative_balance_change": errors,
        "closed_system": config.simulation.closed_system,
        "nonnegative": bool(np.min(final) >= -1e-10),
        "ammonia_free": bool(np.max(final[:, S["ammonium_mol_m3"]]) == 0.0),
        "ph_charge_balance_solved": config.simulation.ph_mode == "dynamic",
        "note": "Open-system balance changes include boundary transport and oxygen exchange.",
    }


def simulate_0d(config: Optional[ModelConfig] = None, geochem: Optional[GeochemLookup] = None) -> SimulationResult:
    config = config or ModelConfig()
    config.validate()
    state = _initial_state(config, np.ones(1), "0d")
    initial = state.copy()
    total_s = config.simulation.days * SECONDS_PER_DAY
    output_s = max(config.simulation.output_interval_days * SECONDS_PER_DAY, 1.0)
    targets = np.unique(np.append(np.arange(0.0, total_s, output_s), total_s))
    frames = [_state_frame(state, 0.0, config, {})]
    time_s = 0.0
    max_step = config.simulation.reaction_step_h * 3600.0
    for target in targets[1:]:
        while time_s < target - 1e-8:
            dt = min(max_step, target - time_s)
            state = _reaction_step(state, time_s, dt, config, geochem)
            time_s += dt
        frames.append(_state_frame(state, target / SECONDS_PER_DAY, config, {}))
    return SimulationResult("0d", pd.concat(frames, ignore_index=True), _summary(state, config, "0d"), _diagnostics(initial, state, config, "0d"), config)


def _capsule_profile_1d(config: ModelConfig) -> Tuple[np.ndarray, np.ndarray]:
    trans = config.transport
    x = np.linspace(0.0, trans.crack_length_mm, trans.nx_1d)
    centers = np.linspace(0.08, 0.92, trans.capsule_count_1d) * trans.crack_length_mm
    profile = np.zeros_like(x)
    sigma = max(trans.capsule_spread_mm, 1e-6)
    for center in centers:
        profile += np.exp(-0.5 * ((x - center) / sigma) ** 2)
    # Dose is configured independently; capsule count controls spatial discreteness.
    profile *= 1.0 / max(float(np.mean(profile)), 1e-30)
    return x, profile


def _capsule_profile_2d(config: ModelConfig) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    trans = config.transport
    x = np.linspace(0.0, trans.crack_length_mm, trans.nx_2d)
    y = np.linspace(0.0, trans.crack_width_mm, trans.ny_2d)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    rng = np.random.RandomState(config.simulation.random_seed)
    centers_x = rng.uniform(0.05, 0.95, trans.capsule_count_2d) * trans.crack_length_mm
    centers_y = rng.choice([0.10, 0.90], trans.capsule_count_2d) * trans.crack_width_mm
    profile = np.zeros_like(xx)
    sigma_x = max(trans.capsule_spread_mm, 1e-6)
    sigma_y = max(trans.crack_width_mm / 5.0, 1e-6)
    for center_x, center_y in zip(centers_x, centers_y):
        profile += np.exp(-0.5 * ((xx - center_x) / sigma_x) ** 2 - 0.5 * ((yy - center_y) / sigma_y) ** 2)
    # Exact mean normalization keeps dose fixed while count changes discreteness.
    profile *= 1.0 / max(float(np.mean(profile)), 1e-30)
    return xx, yy, profile


def _diffusion_operator_1d(diffusivity: np.ndarray, dx: float, velocity: float) -> Tuple[sparse.csr_matrix, np.ndarray]:
    n = diffusivity.size
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    boundary = np.zeros(n)
    for i in range(n):
        diagonal = 0.0
        if i == 0:
            coefficient = 2.0 * diffusivity[i] / dx**2
            diagonal -= coefficient
            boundary[i] += coefficient
        else:
            face = 2.0 * diffusivity[i] * diffusivity[i - 1] / max(diffusivity[i] + diffusivity[i - 1], 1e-30)
            coefficient = face / dx**2
            rows.append(i); cols.append(i - 1); data.append(coefficient)
            diagonal -= coefficient
        if i < n - 1:
            face = 2.0 * diffusivity[i] * diffusivity[i + 1] / max(diffusivity[i] + diffusivity[i + 1], 1e-30)
            coefficient = face / dx**2
            rows.append(i); cols.append(i + 1); data.append(coefficient)
            diagonal -= coefficient
        if velocity > 0.0:
            diagonal -= velocity / dx
            if i == 0:
                boundary[i] += velocity / dx
            else:
                rows.append(i); cols.append(i - 1); data.append(velocity / dx)
        rows.append(i); cols.append(i); data.append(diagonal)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n)), boundary


def _diffusion_operator_2d(diffusivity: np.ndarray, dx: float, dy: float, velocity: float) -> Tuple[sparse.csr_matrix, np.ndarray]:
    ny, nx = diffusivity.shape
    n = nx * ny
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    boundary = np.zeros(n)
    for j in range(ny):
        for i in range(nx):
            index = j * nx + i
            diagonal = 0.0
            for di, dj, distance in ((-1, 0, dx), (1, 0, dx), (0, -1, dy), (0, 1, dy)):
                ni, nj = i + di, j + dj
                if ni < 0:
                    coefficient = 2.0 * diffusivity[j, i] / dx**2
                    diagonal -= coefficient
                    boundary[index] += coefficient
                elif ni >= nx or nj < 0 or nj >= ny:
                    continue
                else:
                    other = diffusivity[nj, ni]
                    face = 2.0 * diffusivity[j, i] * other / max(diffusivity[j, i] + other, 1e-30)
                    coefficient = face / distance**2
                    rows.append(index); cols.append(nj * nx + ni); data.append(coefficient)
                    diagonal -= coefficient
            if velocity > 0.0:
                diagonal -= velocity / dx
                if i == 0:
                    boundary[index] += velocity / dx
                else:
                    rows.append(index); cols.append(index - 1); data.append(velocity / dx)
            rows.append(index); cols.append(index); data.append(diagonal)
    return sparse.csr_matrix((data, (rows, cols)), shape=(n, n)), boundary


def _transport_step(
    state: np.ndarray,
    time_s: float,
    dt_s: float,
    config: ModelConfig,
    shape: Tuple[int, ...],
) -> np.ndarray:
    _, wet_factor = _environment(time_s + 0.5 * dt_s, config)
    level = "2d" if len(shape) == 2 else "1d"
    cell_volume, wall_area = _uniform_cell_geometry(config, level, state.shape[0])
    metrics = repair_metrics(state, config, cell_volume, wall_area)
    obstruction = np.maximum(1.0 - metrics["crack_closure_ratio"], 1e-3) ** config.transport.tortuosity_exponent
    base_diffusivity = {
        "lactate": config.transport.diffusivity_lactate_m2_s,
        "oxygen": config.transport.diffusivity_oxygen_m2_s,
        "calcium": config.transport.diffusivity_calcium_m2_s,
        "carbon": config.transport.diffusivity_carbon_m2_s,
    }
    boundary_values = {
        "lactate": 0.0,
        "oxygen": config.environment.oxygen_boundary_mol_m3,
        "calcium": 0.0,
        "carbon": config.environment.inorganic_carbon_boundary_mol_m3,
    }
    identity = sparse.eye(state.shape[0], format="csr")
    updated = state.copy()
    for index, species in DISSOLVED.items():
        local_d = base_diffusivity[species] * wet_factor * obstruction
        if len(shape) == 1:
            dx = config.transport.crack_length_mm / 1000.0 / (shape[0] - 1)
            operator, boundary = _diffusion_operator_1d(local_d, dx, config.transport.advective_velocity_m_s)
        else:
            ny, nx = shape
            dx = config.transport.crack_length_mm / 1000.0 / (nx - 1)
            dy = config.transport.crack_width_mm / 1000.0 / (ny - 1)
            operator, boundary = _diffusion_operator_2d(
                local_d.reshape(shape), dx, dy, config.transport.advective_velocity_m_s
            )
        rhs = state[:, index] + dt_s * boundary * boundary_values[species]
        updated[:, index] = np.maximum(spsolve(identity - dt_s * operator, rhs), 0.0)
    updated[:, S["hydrated_carbon_mol_m3"]] = np.minimum(
        updated[:, S["hydrated_carbon_mol_m3"]], updated[:, S["inorganic_carbon_mol_m3"]]
    )
    return updated


def _spatial_simulation(
    level: str,
    config: ModelConfig,
    geochem: Optional[GeochemLookup],
) -> SimulationResult:
    if level == "1d":
        x, profile = _capsule_profile_1d(config)
        shape = (x.size,)
        coordinates = {"x_mm": x}
    else:
        xx, yy, profile_2d = _capsule_profile_2d(config)
        profile = profile_2d.ravel()
        shape = profile_2d.shape
        coordinates = {"x_mm": xx.ravel(), "y_mm": yy.ravel()}
    state = _initial_state(config, profile.ravel(), level)
    initial = state.copy()
    frames = [_state_frame(state, 0.0, config, coordinates)]
    total_s = config.simulation.days * SECONDS_PER_DAY
    snapshot_days = sorted(set([0.0, min(7.0, config.simulation.days), min(14.0, config.simulation.days), min(21.0, config.simulation.days), config.simulation.days]))
    targets = [day * SECONDS_PER_DAY for day in snapshot_days]
    time_s = 0.0
    max_step = config.simulation.reaction_step_h * 3600.0
    for target in targets[1:]:
        while time_s < target - 1e-8:
            dt = min(max_step, target - time_s)
            state = _reaction_step(state, time_s, dt, config, geochem)
            if not config.simulation.closed_system:
                state = _transport_step(state, time_s, dt, config, shape)
            time_s += dt
        frames.append(_state_frame(state, target / SECONDS_PER_DAY, config, coordinates))
    diagnostics = _diagnostics(initial, state, config, level)
    diagnostics["grid_shape"] = list(shape)
    diagnostics["true_spatial_solver"] = True
    return SimulationResult(level, pd.concat(frames, ignore_index=True), _summary(state, config, level), diagnostics, config)


def simulate_1d(config: Optional[ModelConfig] = None, geochem: Optional[GeochemLookup] = None) -> SimulationResult:
    config = config or ModelConfig()
    config.validate()
    return _spatial_simulation("1d", config, geochem)


def simulate_2d(config: Optional[ModelConfig] = None, geochem: Optional[GeochemLookup] = None) -> SimulationResult:
    config = config or ModelConfig()
    config.validate()
    return _spatial_simulation("2d", config, geochem)

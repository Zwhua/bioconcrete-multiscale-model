"""Configuration and parameter provenance for the BioConcrete model."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Type, TypeVar


@dataclass
class SimulationConfig:
    days: float = 28.0
    output_interval_days: float = 1.0
    reaction_step_h: float = 12.0
    carbonate_mode: str = "equilibrium"
    ph_mode: str = "dynamic"
    closed_system: bool = False
    random_seed: int = 2026


@dataclass
class EnvironmentConfig:
    water_activity_wet: float = 1.0
    water_activity_dry: float = 0.90
    oxygen_boundary_mol_m3: float = 0.25
    oxygen_initial_mol_m3: float = 0.02
    inorganic_carbon_boundary_mol_m3: float = 0.015
    inorganic_carbon_initial_mol_m3: float = 0.0
    ph: float = 11.5
    temperature_c: float = 30.0
    exposure: str = "intermittent"
    wet_hours_per_day: float = 12.0
    oxygen_transfer_s: float = 2.0e-6
    initial_alkalinity_mol_m3: float = 35.0
    ph_minimum: float = 6.0
    ph_maximum: float = 13.5


@dataclass
class MicrobialKineticsConfig:
    """Anonymous population-scale activity without construction details."""

    agent_dosage_multiplier: float = 1.0
    dosage_basis: str = "fixed_total_inventory"
    reference_inventory_volume_m3: float = 6.0e-7
    capsule_calcium_lactate_mol_m3: float = 8000.0
    capsule_csh_volume_fraction: float = 0.035
    spore_density_rel: float = 1.0
    active_density_rel: float = 0.0
    capsule_release_s: float = 1.3e-6
    csh_release_s: float = 1.3e-6
    germination_s: float = 7.0e-6
    maximum_growth_s: float = 7.0e-6
    decay_s: float = 1.5e-7
    alkaline_decay_s: float = 4.0e-7
    encapsulation_decay_m3_mol: float = 1.0e-5
    qmax_lactate_mol_m3_s: float = 2.5e-3
    effective_kcat_s: float = 2.5e-3
    effective_km_mol_m3: float = 350.0
    active_unit_concentration: float = 1.0
    activity_multiplier: float = 1.0
    k_lactate_mol_m3: float = 350.0
    k_oxygen_mol_m3: float = 0.035
    biomass_carbon_fraction: float = 0.08
    biomass_yield_rel_m3_mol: float = 1.0e-4
    aw_threshold: float = 0.95
    oxygen_threshold_mol_m3: float = 0.06
    ph_optimum: float = 10.5
    ph_width: float = 2.0
    temperature_optimum_c: float = 30.0
    temperature_width_c: float = 12.0
    response_delay_h: float = 4.0
    basal_leak_fraction: float = 0.01
    activation_duration_h: float = 4.0
    signal_relaxation_h: float = 1.0
    oxygen_rise_threshold_mol_m3_h: float = 0.005
    ph_drop_threshold_h: float = 0.05
    gate_logic: str = "AND"


@dataclass
class ChemistryConfig:
    portlandite_mol_m3: float = 5000.0
    calcite_rate_mol_m3_s: float = 1.0e-4
    calcite_surface_area_rel: float = 1.0
    calcite_reaction_order: float = 1.0
    portlandite_dissolution_s: float = 2.0e-7
    portlandite_carbon_half_mol_m3: float = 50.0
    calcite_ksp: float = 3.31e-9
    ca_hydration_s: float = 2.0e-3
    ca_dehydration_s: float = 5.0e-4
    equilibrium_relaxation_s: float = 2.0e-2
    calcite_molar_mass_kg_mol: float = 0.1000869
    calcite_density_kg_m3: float = 2710.0
    csh_density_kg_m3: float = 2400.0
    wall_deposition_fraction: float = 0.75


@dataclass
class TransportConfig:
    crack_length_mm: float = 100.0
    crack_width_mm: float = 0.30
    crack_depth_mm: float = 20.0
    out_of_plane_thickness_mm: float = 1.0
    nx_1d: int = 51
    nx_2d: int = 15
    ny_2d: int = 5
    capsule_count_1d: int = 8
    capsule_count_2d: int = 8
    capsule_spread_mm: float = 2.0
    diffusivity_lactate_m2_s: float = 7.0e-10
    diffusivity_oxygen_m2_s: float = 1.9e-9
    diffusivity_calcium_m2_s: float = 0.79e-9
    diffusivity_carbon_m2_s: float = 1.9e-9
    dry_diffusivity_factor: float = 0.02
    porosity_initial: float = 0.18
    porosity_minimum: float = 0.02
    tortuosity_exponent: float = 1.5
    advective_velocity_m_s: float = 0.0


@dataclass
class Geometry3DConfig:
    mode: str = "rectangular"
    topology: str = "blind_crack"
    nx: int = 51
    ny: int = 5
    nz: int = 21
    capsule_count: int = 24
    capsule_spread_x_mm: float = 2.0
    capsule_spread_y_mm: float = 0.06
    capsule_spread_z_mm: float = 1.0
    capsule_depth_mode: str = "uniform"
    aperture_field_path: Optional[str] = None


@dataclass
class Boundary3DConfig:
    x_min: str = "exposed"
    x_max: str = "no_flux"
    y_min: str = "crack_wall"
    y_max: str = "crack_wall"
    z_min: str = "no_flux"
    z_max: str = "no_flux"
    oxygen_supply_mode: str = "boundary_robin"
    oxygen_mass_transfer_m_s: float = 2.0e-6
    carbon_mass_transfer_m_s: float = 2.0e-6


@dataclass
class Solver3DConfig:
    splitting_scheme: str = "strang"
    linear_solver: str = "auto"
    reaction_batch_cells: int = 64
    relative_tolerance: float = 1.0e-6
    absolute_tolerance: float = 1.0e-10
    maximum_linear_iterations: int = 500
    checkpoint_interval_steps: int = 24
    memory_limit_gb: float = 8.0
    reaction_workers: int = 1
    reaction_parallel_backend: str = "serial"


@dataclass
class Output3DConfig:
    storage_format: str = "zarr"
    storage_dtype: str = "float64"
    save_full_state: bool = True
    save_every_days: float = 1.0
    compression: str = "zstd"


@dataclass
class ModelConfig:
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    kinetics: MicrobialKineticsConfig = field(default_factory=MicrobialKineticsConfig)
    chemistry: ChemistryConfig = field(default_factory=ChemistryConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)
    geometry_3d: Geometry3DConfig = field(default_factory=Geometry3DConfig)
    boundary_3d: Boundary3DConfig = field(default_factory=Boundary3DConfig)
    solver_3d: Solver3DConfig = field(default_factory=Solver3DConfig)
    output_3d: Output3DConfig = field(default_factory=Output3DConfig)

    def validate(self) -> None:
        if self.simulation.days <= 0 or self.simulation.reaction_step_h <= 0:
            raise ValueError("Simulation duration and time step must be positive")
        if self.simulation.carbonate_mode not in {"equilibrium", "kinetic"}:
            raise ValueError("carbonate_mode must be 'equilibrium' or 'kinetic'")
        if self.simulation.ph_mode not in {"dynamic", "fixed"}:
            raise ValueError("ph_mode must be 'dynamic' or 'fixed'")
        if self.environment.exposure not in {"intermittent", "continuous", "dry"}:
            raise ValueError("exposure must be intermittent, continuous, or dry")
        if not 0 <= self.kinetics.biomass_carbon_fraction < 1:
            raise ValueError("biomass_carbon_fraction must be in [0, 1)")
        if not 0 <= self.kinetics.basal_leak_fraction <= 1:
            raise ValueError("basal_leak_fraction must be in [0, 1]")
        if self.kinetics.agent_dosage_multiplier < 0:
            raise ValueError("agent_dosage_multiplier must be nonnegative")
        if self.kinetics.dosage_basis not in {"fixed_total_inventory", "fixed_concentration"}:
            raise ValueError("dosage_basis must be fixed_total_inventory or fixed_concentration")
        if self.kinetics.reference_inventory_volume_m3 <= 0:
            raise ValueError("reference_inventory_volume_m3 must be positive")
        if self.kinetics.gate_logic not in {"AND", "OR", "static_suitability"}:
            raise ValueError("gate_logic must be AND, OR, or static_suitability")
        if not 0 <= self.chemistry.wall_deposition_fraction <= 1:
            raise ValueError("wall_deposition_fraction must be in [0, 1]")
        if min(self.transport.nx_1d, self.transport.nx_2d, self.transport.ny_2d) < 3:
            raise ValueError("Each spatial grid dimension must contain at least three cells")
        if self.transport.crack_width_mm <= 0 or self.transport.crack_length_mm <= 0:
            raise ValueError("Crack dimensions must be positive")
        if self.geometry_3d.mode not in {"rectangular", "aperture_field", "voxel_ct"}:
            raise ValueError("Unsupported 3D geometry mode")
        if self.geometry_3d.topology not in {"blind_crack", "through_crack"}:
            raise ValueError("Unsupported crack topology")
        if min(self.geometry_3d.nx, self.geometry_3d.ny, self.geometry_3d.nz) < 2:
            raise ValueError("Each 3D grid dimension must contain at least two cells")
        if self.geometry_3d.capsule_depth_mode not in {"uniform", "surface", "layered"}:
            raise ValueError("capsule_depth_mode must be uniform, surface, or layered")
        if self.boundary_3d.oxygen_supply_mode not in {"legacy_volumetric", "boundary_robin"}:
            raise ValueError("Unsupported 3D oxygen supply mode")
        if self.solver_3d.splitting_scheme not in {"lie", "strang"}:
            raise ValueError("splitting_scheme must be lie or strang")
        if self.solver_3d.linear_solver not in {"direct", "cg", "gmres", "bicgstab", "auto"}:
            raise ValueError("Unsupported 3D linear solver")
        if self.solver_3d.reaction_batch_cells <= 0 or self.solver_3d.memory_limit_gb <= 0:
            raise ValueError("3D solver batch and memory limits must be positive")
        if self.solver_3d.reaction_workers <= 0:
            raise ValueError("reaction_workers must be positive")
        if self.solver_3d.reaction_parallel_backend not in {"serial", "thread", "process"}:
            raise ValueError("reaction_parallel_backend must be serial, thread, or process")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Optional[Path] = None) -> "ModelConfig":
        if path is None:
            config = cls()
        else:
            raw = json.loads(path.read_text(encoding="utf-8"))
            config = cls(
                simulation=_load_dataclass(SimulationConfig, raw.get("simulation", {})),
                environment=_load_dataclass(EnvironmentConfig, raw.get("environment", {})),
                kinetics=_load_dataclass(MicrobialKineticsConfig, raw.get("kinetics", {})),
                chemistry=_load_dataclass(ChemistryConfig, raw.get("chemistry", {})),
                transport=_load_dataclass(TransportConfig, raw.get("transport", {})),
                geometry_3d=_load_dataclass(Geometry3DConfig, raw.get("geometry_3d", {})),
                boundary_3d=_load_dataclass(Boundary3DConfig, raw.get("boundary_3d", {})),
                solver_3d=_load_dataclass(Solver3DConfig, raw.get("solver_3d", {})),
                output_3d=_load_dataclass(Output3DConfig, raw.get("output_3d", {})),
            )
        config.validate()
        return config


T = TypeVar("T")


def _load_dataclass(kind: Type[T], values: Mapping[str, Any]) -> T:
    allowed = {item.name for item in fields(kind)}
    unknown = set(values) - allowed
    if unknown:
        raise ValueError("Unknown {} settings: {}".format(kind.__name__, sorted(unknown)))
    return kind(**dict(values))


PARAMETER_PROVENANCE = {
    "kinetics.agent_dosage_multiplier": (
        "scenario variable", "dimensionless complete repair-agent dose multiplier", 0.0, 3.0
    ),
    "kinetics.qmax_lactate_mol_m3_s": ("model prior", "aggregate substrate utilization rate; requires calibration", 5.0e-5, 1.0e-2),
    "kinetics.k_lactate_mol_m3": ("database aggregate", "anonymous summary of public kinetic records", 1.0, 5000.0),
    "kinetics.k_oxygen_mol_m3": ("literature prior", "aerobic Monod half-saturation", 0.005, 0.20),
    "kinetics.decay_s": ("literature prior", "protected population-scale decay", 1.0e-8, 2.0e-6),
    "kinetics.capsule_release_s": ("project hypothesis", "effective release rate", 1.0e-7, 3.0e-5),
    "chemistry.calcite_rate_mol_m3_s": ("literature prior", "saturation-index precipitation law", 1.0e-7, 1.0e-2),
    "chemistry.calcite_ksp": ("thermodynamic database", "PHREEQC calcite at ambient temperature", 2.5e-9, 5.0e-9),
    "transport.diffusivity_oxygen_m2_s": ("literature prior", "aqueous oxygen/CO2 scale", 5.0e-10, 3.0e-9),
    "transport.crack_width_mm": ("preregistered scenario", "fixed prospective crack range", 0.05, 0.50),
    "chemistry.portlandite_mol_m3": ("project hypothesis", "effective crack-adjacent reservoir", 100.0, 15000.0),
    "kinetics.effective_kcat_s": ("database aggregate", "anonymous effective catalytic prior", 5.0e-5, 1.0e-2),
    "kinetics.response_delay_h": ("scenario prior", "anonymous response delay", 0.0, 24.0),
    "kinetics.basal_leak_fraction": ("scenario prior", "anonymous basal activity", 0.0, 0.10),
    "chemistry.wall_deposition_fraction": ("literature prior", "fraction of solid deposited at crack walls; public calibration pending", 0.05, 1.0),
    "kinetics.activity_multiplier": ("scenario prior", "anonymous effective activity; public calibration pending", 0.5, 5.0),
}

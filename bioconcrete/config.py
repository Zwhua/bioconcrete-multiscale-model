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
    closed_system: bool = False
    random_seed: int = 2026


@dataclass
class EnvironmentConfig:
    water_activity_wet: float = 1.0
    water_activity_dry: float = 0.90
    oxygen_boundary_mol_m3: float = 0.25
    oxygen_initial_mol_m3: float = 0.02
    inorganic_carbon_boundary_mol_m3: float = 0.015
    ph: float = 11.5
    temperature_c: float = 30.0
    exposure: str = "intermittent"
    wet_hours_per_day: float = 12.0
    oxygen_transfer_s: float = 2.0e-6


@dataclass
class MicrobialKineticsConfig:
    """Population-scale kinetics without strain, sequence, or circuit details."""

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


@dataclass
class TransportConfig:
    crack_length_mm: float = 100.0
    crack_width_mm: float = 0.30
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
class ModelConfig:
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    kinetics: MicrobialKineticsConfig = field(default_factory=MicrobialKineticsConfig)
    chemistry: ChemistryConfig = field(default_factory=ChemistryConfig)
    transport: TransportConfig = field(default_factory=TransportConfig)

    def validate(self) -> None:
        if self.simulation.days <= 0 or self.simulation.reaction_step_h <= 0:
            raise ValueError("Simulation duration and time step must be positive")
        if self.simulation.carbonate_mode not in {"equilibrium", "kinetic"}:
            raise ValueError("carbonate_mode must be 'equilibrium' or 'kinetic'")
        if self.environment.exposure not in {"intermittent", "continuous", "dry"}:
            raise ValueError("exposure must be intermittent, continuous, or dry")
        if not 0 <= self.kinetics.biomass_carbon_fraction < 1:
            raise ValueError("biomass_carbon_fraction must be in [0, 1)")
        if min(self.transport.nx_1d, self.transport.nx_2d, self.transport.ny_2d) < 3:
            raise ValueError("Each spatial grid dimension must contain at least three cells")
        if self.transport.crack_width_mm <= 0 or self.transport.crack_length_mm <= 0:
            raise ValueError("Crack dimensions must be positive")

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
    "kinetics.qmax_lactate_mol_m3_s": ("model prior", "aggregate substrate utilization rate; requires calibration", 5.0e-5, 1.0e-2),
    "kinetics.k_lactate_mol_m3": ("database aggregate", "public kinetic records summarized without sequence information", 1.0, 5000.0),
    "kinetics.k_oxygen_mol_m3": ("literature prior", "aerobic Monod half-saturation", 0.005, 0.20),
    "kinetics.decay_s": ("literature prior", "protected population-scale decay", 1.0e-8, 2.0e-6),
    "kinetics.capsule_release_s": ("project hypothesis", "effective release rate", 1.0e-7, 3.0e-5),
    "chemistry.calcite_rate_mol_m3_s": ("literature prior", "saturation-index precipitation law", 1.0e-7, 1.0e-2),
    "chemistry.calcite_ksp": ("thermodynamic database", "PHREEQC calcite at ambient temperature", 2.5e-9, 5.0e-9),
    "transport.diffusivity_oxygen_m2_s": ("literature prior", "aqueous oxygen/CO2 scale", 5.0e-10, 3.0e-9),
    "transport.crack_width_mm": ("project experiment", "planned crack range 0.1-0.5 mm", 0.05, 0.50),
    "chemistry.portlandite_mol_m3": ("project hypothesis", "effective crack-adjacent reservoir", 100.0, 15000.0),
}

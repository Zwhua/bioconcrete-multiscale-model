"""Canonical local-state schema shared by every model dimension."""

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

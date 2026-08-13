"""Typed, species-specific boundary conditions for the 3D finite-volume model."""

from dataclasses import dataclass
from typing import Dict


BOUNDARY_KINDS = {"no_flux", "dirichlet", "robin", "inlet", "outlet", "crack_wall"}
FACES = {"x_min", "x_max", "y_min", "y_max", "z_min", "z_max"}


@dataclass(frozen=True)
class BoundaryCondition3D:
    kind: str = "no_flux"
    value: float = 0.0
    mass_transfer_m_s: float = 0.0

    def __post_init__(self):
        if self.kind not in BOUNDARY_KINDS:
            raise ValueError("Unsupported boundary kind: {}".format(self.kind))
        if self.mass_transfer_m_s < 0:
            raise ValueError("mass_transfer_m_s must be nonnegative")


def no_flux_boundaries() -> Dict[str, BoundaryCondition3D]:
    return {face: BoundaryCondition3D() for face in FACES}


def validate_boundaries(boundaries: Dict[str, BoundaryCondition3D], topology: str) -> None:
    unknown = set(boundaries) - FACES
    if unknown:
        raise ValueError("Unknown 3D faces: {}".format(sorted(unknown)))
    merged = no_flux_boundaries()
    merged.update(boundaries)
    if topology == "blind_crack" and merged["x_max"].kind in {"outlet", "inlet"}:
        raise ValueError("A blind crack cannot have an x_max inlet or outlet")
    if topology == "through_crack":
        if merged["x_min"].kind != "inlet" or merged["x_max"].kind != "outlet":
            raise ValueError("A through crack requires explicit x_min inlet and x_max outlet")

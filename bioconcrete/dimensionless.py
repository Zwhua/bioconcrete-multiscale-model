"""Dimensionless screening for choosing the 2.5D or full 3D model.

The calculations are deliberately independent of the solver so they can be run
before allocating a three-dimensional field.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Dict, Optional

from .config import ModelConfig


@dataclass(frozen=True)
class DimensionlessSummary:
    aspect_length_to_aperture: float
    aspect_depth_to_aperture: float
    diffusion_time_aperture_s: float
    diffusion_time_length_s: float
    diffusion_time_depth_s: float
    advection_time_length_s: Optional[float]
    reaction_time_s: float
    wet_period_s: float
    peclet_length: float
    damkohler_length: float
    reynolds_aperture: float
    reaction_to_wet_period: float
    aperture_to_reaction_ratio: float
    aperture_to_wet_period_ratio: float
    aperture_to_inplane_transport_ratio: float
    two_point_five_d_applicable: bool
    applicability_reason: str
    evidence_label: str = "uncalibrated 3D model output; not experimental data"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def analyze_dimensionless(
    config: ModelConfig,
    *,
    diffusivity_m2_s: float = None,
    reaction_rate_s: float = None,
    density_kg_m3: float = 1000.0,
    dynamic_viscosity_pa_s: float = 1.0e-3,
    separation_factor: float = 0.1,
) -> DimensionlessSummary:
    """Return scale analysis using characteristic crack dimensions.

    ``separation_factor`` makes "far smaller" operational: all aperture mixing
    ratios must be below it for the 2.5D screen to pass.
    """
    if diffusivity_m2_s is None:
        diffusivity_m2_s = config.transport.diffusivity_oxygen_m2_s
    if reaction_rate_s is None:
        reaction_rate_s = config.kinetics.maximum_growth_s
    values = (diffusivity_m2_s, reaction_rate_s, density_kg_m3,
              dynamic_viscosity_pa_s, separation_factor)
    if any(value <= 0 for value in values):
        raise ValueError("Dimensionless inputs must be positive")

    length = config.transport.crack_length_mm * 1.0e-3
    aperture = config.transport.crack_width_mm * 1.0e-3
    depth = config.transport.crack_depth_mm * 1.0e-3
    velocity = abs(config.transport.advective_velocity_m_s)

    tau_y = aperture ** 2 / diffusivity_m2_s
    tau_x = length ** 2 / diffusivity_m2_s
    tau_z = depth ** 2 / diffusivity_m2_s
    tau_reaction = 1.0 / reaction_rate_s
    wet_period = config.environment.wet_hours_per_day * 3600.0
    advection_time = length / velocity if velocity > 0 else None
    inplane_time = min(tau_x, tau_z, advection_time if advection_time is not None else tau_x)
    ratios = (tau_y / tau_reaction, tau_y / wet_period, tau_y / inplane_time)
    applicable = all(value < separation_factor for value in ratios)
    reason = (
        "aperture mixing is separated from reaction, wetting, and in-plane transport"
        if applicable else
        "aperture mixing is not sufficiently faster than every competing timescale; use full 3D"
    )
    return DimensionlessSummary(
        aspect_length_to_aperture=length / aperture,
        aspect_depth_to_aperture=depth / aperture,
        diffusion_time_aperture_s=tau_y,
        diffusion_time_length_s=tau_x,
        diffusion_time_depth_s=tau_z,
        advection_time_length_s=advection_time,
        reaction_time_s=tau_reaction,
        wet_period_s=wet_period,
        peclet_length=velocity * length / diffusivity_m2_s,
        damkohler_length=reaction_rate_s * length ** 2 / diffusivity_m2_s,
        reynolds_aperture=density_kg_m3 * velocity * aperture / dynamic_viscosity_pa_s,
        reaction_to_wet_period=tau_reaction / wet_period,
        aperture_to_reaction_ratio=ratios[0],
        aperture_to_wet_period_ratio=ratios[1],
        aperture_to_inplane_transport_ratio=ratios[2],
        two_point_five_d_applicable=applicable,
        applicability_reason=reason,
    )


def write_dimensionless_summary(
    summary: DimensionlessSummary,
    output: Path = Path("model_runs/v0.6.0/dimensionless/dimensionless_summary.json"),
) -> Path:
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(summary.to_dict(), indent=2), encoding="utf-8")
    temporary.replace(output)
    return output

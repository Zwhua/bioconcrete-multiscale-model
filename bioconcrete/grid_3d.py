"""Cell-centred structured 3D crack geometry in SI units."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Dict, Tuple

import numpy as np

from .config import ModelConfig


@dataclass(frozen=True)
class StructuredGrid3D:
    x_m: np.ndarray
    y_m: np.ndarray
    z_m: np.ndarray
    cell_volume_m3: np.ndarray
    face_area_x_m2: np.ndarray
    face_area_y_m2: np.ndarray
    face_area_z_m2: np.ndarray
    fluid_mask: np.ndarray

    @property
    def shape(self) -> Tuple[int, int, int]:
        return (self.z_m.size, self.y_m.size, self.x_m.size)

    @property
    def size(self) -> int:
        return int(np.prod(self.shape))

    def flatten(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if values.shape != self.shape:
            raise ValueError("field shape must be {}".format(self.shape))
        return values.ravel(order="C")

    def reshape(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values)
        if values.size != self.size:
            raise ValueError("flat field size must be {}".format(self.size))
        return values.reshape(self.shape, order="C")

    def cell_coordinates(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        return np.meshgrid(self.z_m, self.y_m, self.x_m, indexing="ij")

    def boundary_indices(self) -> Dict[str, np.ndarray]:
        ids = np.arange(self.size).reshape(self.shape)
        return {"x_min": ids[:, :, 0].ravel(), "x_max": ids[:, :, -1].ravel(),
                "y_min": ids[:, 0, :].ravel(), "y_max": ids[:, -1, :].ravel(),
                "z_min": ids[0, :, :].ravel(), "z_max": ids[-1, :, :].ravel()}

    @property
    def geometry_hash(self) -> str:
        digest = hashlib.sha256()
        for value in (self.x_m, self.y_m, self.z_m, self.cell_volume_m3,
                      self.face_area_x_m2, self.face_area_y_m2,
                      self.face_area_z_m2, self.fluid_mask):
            digest.update(np.ascontiguousarray(value).tobytes())
        digest.update(b"axis_order=z,y,x;flatten_order=C")
        return digest.hexdigest()


def rectangular_grid_3d(config: ModelConfig) -> StructuredGrid3D:
    geometry = config.geometry_3d
    if geometry.mode != "rectangular":
        if not geometry.aperture_field_path:
            raise ValueError("{} geometry requires an explicit data path".format(geometry.mode))
        raise NotImplementedError("validated {} reader is reserved for a later version".format(geometry.mode))
    nx, ny, nz = geometry.nx, geometry.ny, geometry.nz
    lx = config.transport.crack_length_mm * 1e-3
    ly = config.transport.crack_width_mm * 1e-3
    lz = config.transport.crack_depth_mm * 1e-3
    dx, dy, dz = lx / nx, ly / ny, lz / nz
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    z = (np.arange(nz) + 0.5) * dz
    shape = (nz, ny, nx)
    return StructuredGrid3D(
        x, y, z, np.full(shape, dx * dy * dz),
        np.full(shape, dy * dz), np.full(shape, dx * dz),
        np.full(shape, dx * dy), np.ones(shape, dtype=bool),
    )


def capsule_profile_3d(config: ModelConfig, grid: StructuredGrid3D) -> Tuple[np.ndarray, np.ndarray]:
    """Return volume-normalised profile and physical capsule centres `(x,y,z)`."""
    g = config.geometry_3d
    rng = np.random.RandomState(config.simulation.random_seed)
    count = g.capsule_count
    lx = (grid.x_m[1] - grid.x_m[0]) * grid.x_m.size
    ly = (grid.y_m[1] - grid.y_m[0]) * grid.y_m.size
    lz = (grid.z_m[1] - grid.z_m[0]) * grid.z_m.size
    cx = rng.uniform(0.05 * lx, 0.95 * lx, count)
    cy = rng.choice([0.1, 0.9], count) * ly
    if g.capsule_depth_mode == "surface":
        cz = rng.uniform(0.0, 0.2 * lz, count)
    elif g.capsule_depth_mode == "layered":
        cz = rng.choice([0.2, 0.5, 0.8], count) * lz
    else:
        cz = rng.uniform(0.05 * lz, 0.95 * lz, count)
    zz, yy, xx = grid.cell_coordinates()
    profile = np.zeros(grid.shape)
    sx, sy, sz = (max(g.capsule_spread_x_mm * 1e-3, 1e-12),
                  max(g.capsule_spread_y_mm * 1e-3, 1e-12),
                  max(g.capsule_spread_z_mm * 1e-3, 1e-12))
    for x0, y0, z0 in zip(cx, cy, cz):
        profile += np.exp(-0.5 * (((xx-x0)/sx)**2 + ((yy-y0)/sy)**2 + ((zz-z0)/sz)**2))
    weighted_mean = np.sum(profile * grid.cell_volume_m3) / np.sum(grid.cell_volume_m3)
    if weighted_mean <= 0:
        raise ValueError("capsule profile has zero support on this grid")
    return profile / weighted_mean, np.column_stack((cx, cy, cz))

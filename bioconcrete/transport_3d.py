"""Conservative cell-centred seven-point finite-volume transport in 3D."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Union

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import bicgstab, cg, gmres, spsolve

from .boundary_3d import BoundaryCondition3D, no_flux_boundaries
from .grid_3d import StructuredGrid3D


@dataclass(frozen=True)
class LinearSolveDiagnostics:
    solver: str
    converged: bool
    iterations: int
    residual_norm: float
    reason: str


@dataclass(frozen=True)
class TransportDiagnostics3D:
    linear: LinearSolveDiagnostics
    boundary_rate_before: float
    boundary_rate_after: float
    boundary_rates_after: Dict[str, float]


def _spacing(values: np.ndarray) -> float:
    if values.size < 2:
        raise ValueError("3D axes require at least two cells")
    difference = np.diff(values)
    if not np.allclose(difference, difference[0]):
        raise NotImplementedError("The first 3D solver supports uniform structured axes")
    return float(difference[0])


def build_transport_operator_3d(
    grid: StructuredGrid3D,
    diffusivity_m2_s: Union[float, np.ndarray],
    boundaries: Optional[Dict[str, BoundaryCondition3D]] = None,
    velocity_x_m_s: float = 0.0,
) -> Tuple[sparse.csr_matrix, np.ndarray]:
    """Build ``dC/dt = L C + source`` using physical faces and volumes."""
    bc = no_flux_boundaries()
    if boundaries:
        bc.update(boundaries)
    shape = grid.shape
    d = np.broadcast_to(np.asarray(diffusivity_m2_s, dtype=float), shape)
    if np.any(d < 0):
        raise ValueError("diffusivity must be nonnegative")
    dx, dy, dz = _spacing(grid.x_m), _spacing(grid.y_m), _spacing(grid.z_m)
    rows, cols, data = [], [], []
    source = np.zeros(grid.size)
    ids = np.arange(grid.size).reshape(shape)
    specs = ((2, dx, grid.face_area_x_m2), (1, dy, grid.face_area_y_m2),
             (0, dz, grid.face_area_z_m2))
    for axis, distance, areas in specs:
        left = [slice(None)] * 3; right = [slice(None)] * 3
        left[axis] = slice(0, -1); right[axis] = slice(1, None)
        left, right = tuple(left), tuple(right)
        dl, dr = d[left], d[right]
        harmonic = np.where(dl + dr > 0, 2 * dl * dr / np.maximum(dl + dr, 1e-300), 0.0)
        conductance = harmonic * areas[left] / distance
        for p, n, conduct in zip(ids[left].ravel(), ids[right].ravel(), conductance.ravel()):
            if not (grid.fluid_mask.ravel()[p] and grid.fluid_mask.ravel()[n]):
                continue
            cp = conduct / grid.cell_volume_m3.ravel()[p]
            cn = conduct / grid.cell_volume_m3.ravel()[n]
            rows.extend((p, p, n, n)); cols.extend((n, p, p, n)); data.extend((cp, -cp, cn, -cn))

    face_specs = {
        "x_min": (ids[:, :, 0], d[:, :, 0], grid.face_area_x_m2[:, :, 0], dx/2),
        "x_max": (ids[:, :, -1], d[:, :, -1], grid.face_area_x_m2[:, :, -1], dx/2),
        "y_min": (ids[:, 0, :], d[:, 0, :], grid.face_area_y_m2[:, 0, :], dy/2),
        "y_max": (ids[:, -1, :], d[:, -1, :], grid.face_area_y_m2[:, -1, :], dy/2),
        "z_min": (ids[0, :, :], d[0, :, :], grid.face_area_z_m2[0, :, :], dz/2),
        "z_max": (ids[-1, :, :], d[-1, :, :], grid.face_area_z_m2[-1, :, :], dz/2),
    }
    for face, (indices, local_d, areas, half_distance) in face_specs.items():
        condition = bc[face]
        if condition.kind in {"no_flux", "crack_wall", "outlet"}:
            continue
        if condition.kind in {"dirichlet", "inlet"}:
            conductance = local_d * areas / half_distance
        elif condition.kind == "robin":
            conductance = condition.mass_transfer_m_s * areas
        else:
            continue
        for index, conduct in zip(indices.ravel(), conductance.ravel()):
            coefficient = conduct / grid.cell_volume_m3.ravel()[index]
            rows.append(index); cols.append(index); data.append(-coefficient)
            source[index] += coefficient * condition.value

    # Conservative positive-x upwind flux. No advective flux at a no-flux tip.
    if velocity_x_m_s != 0:
        if velocity_x_m_s < 0:
            raise NotImplementedError("The initial solver supports nonnegative x velocity")
        area = grid.face_area_x_m2
        for k in range(shape[0]):
            for j in range(shape[1]):
                for i in range(1, shape[2]):
                    p, n = ids[k,j,i-1], ids[k,j,i]
                    flux = velocity_x_m_s * area[k,j,i-1]
                    rows.extend((p,n)); cols.extend((p,p))
                    data.extend((-flux/grid.cell_volume_m3[k,j,i-1], flux/grid.cell_volume_m3[k,j,i]))
        inlet = bc["x_min"]
        if inlet.kind == "inlet":
            for index, face_area in zip(ids[:,:,0].ravel(), area[:,:,0].ravel()):
                source[index] += velocity_x_m_s * face_area / grid.cell_volume_m3.ravel()[index] * inlet.value
        outlet = bc["x_max"]
        if outlet.kind == "outlet":
            for index, face_area in zip(ids[:,:,-1].ravel(), area[:,:,-1].ravel()):
                rows.append(index); cols.append(index)
                data.append(-velocity_x_m_s * face_area / grid.cell_volume_m3.ravel()[index])
    return sparse.csr_matrix((data, (rows, cols)), shape=(grid.size, grid.size)), source


def _solve(matrix, rhs, method, rtol, atol, maximum_iterations):
    chosen = "direct" if method == "auto" and matrix.shape[0] < 5000 else ("cg" if method == "auto" else method)
    iterations = [0]
    def callback(_): iterations[0] += 1
    if chosen == "direct":
        solution = spsolve(matrix, rhs); info = 0
    else:
        function = {"cg": cg, "gmres": gmres, "bicgstab": bicgstab}[chosen]
        kwargs = {"maxiter": maximum_iterations, "callback": callback}
        try:
            solution, info = function(matrix, rhs, rtol=rtol, atol=atol, **kwargs)
        except TypeError:  # SciPy 1.10 / Python 3.8 compatibility
            solution, info = function(matrix, rhs, tol=rtol, atol=atol, **kwargs)
    residual = float(np.linalg.norm(matrix @ solution - rhs))
    if info != 0 or not np.all(np.isfinite(solution)):
        raise RuntimeError("{} failed (info={}, residual={})".format(chosen, info, residual))
    return solution, LinearSolveDiagnostics(chosen, True, iterations[0], residual, "converged")


def transport_step_3d(
    concentration: np.ndarray, grid: StructuredGrid3D, dt_s: float,
    diffusivity_m2_s: Union[float, np.ndarray],
    boundaries: Optional[Dict[str, BoundaryCondition3D]] = None,
    velocity_x_m_s: float = 0.0, linear_solver: str = "auto",
    relative_tolerance: float = 1e-8, absolute_tolerance: float = 1e-12,
    maximum_iterations: int = 500,
) -> Tuple[np.ndarray, TransportDiagnostics3D]:
    if dt_s <= 0: raise ValueError("dt_s must be positive")
    initial = grid.flatten(np.asarray(concentration, dtype=float))
    operator, source = build_transport_operator_3d(grid, diffusivity_m2_s, boundaries, velocity_x_m_s)
    matrix = sparse.eye(grid.size, format="csr") - dt_s * operator
    solution, diagnostics = _solve(matrix, initial + dt_s * source, linear_solver,
                                   relative_tolerance, absolute_tolerance, maximum_iterations)
    if np.min(solution) < -max(absolute_tolerance * 10, 1e-12):
        raise RuntimeError("transport produced a negative concentration")
    solution = np.maximum(solution, 0.0)
    rate_before = float(np.sum((operator @ initial + source) * grid.cell_volume_m3.ravel()))
    rate_after = float(np.sum((operator @ solution + source) * grid.cell_volume_m3.ravel()))
    face_rates = boundary_flux_rates_3d(
        grid.reshape(solution), grid, diffusivity_m2_s, boundaries, velocity_x_m_s
    )
    return grid.reshape(solution), TransportDiagnostics3D(
        diagnostics, rate_before, rate_after, face_rates
    )


def boundary_flux_rates_3d(
    concentration: np.ndarray,
    grid: StructuredGrid3D,
    diffusivity_m2_s: Union[float, np.ndarray],
    boundaries: Optional[Dict[str, BoundaryCondition3D]] = None,
    velocity_x_m_s: float = 0.0,
) -> Dict[str, float]:
    """Return signed rates into the domain for every boundary face."""
    bc = no_flux_boundaries()
    if boundaries:
        bc.update(boundaries)
    c = np.asarray(concentration, dtype=float)
    d = np.broadcast_to(np.asarray(diffusivity_m2_s, dtype=float), grid.shape)
    dx, dy, dz = _spacing(grid.x_m), _spacing(grid.y_m), _spacing(grid.z_m)
    specifications = {
        "x_min": (c[:, :, 0], d[:, :, 0], grid.face_area_x_m2[:, :, 0], dx / 2),
        "x_max": (c[:, :, -1], d[:, :, -1], grid.face_area_x_m2[:, :, -1], dx / 2),
        "y_min": (c[:, 0, :], d[:, 0, :], grid.face_area_y_m2[:, 0, :], dy / 2),
        "y_max": (c[:, -1, :], d[:, -1, :], grid.face_area_y_m2[:, -1, :], dy / 2),
        "z_min": (c[0, :, :], d[0, :, :], grid.face_area_z_m2[0, :, :], dz / 2),
        "z_max": (c[-1, :, :], d[-1, :, :], grid.face_area_z_m2[-1, :, :], dz / 2),
    }
    rates = {face: 0.0 for face in specifications}
    for face, (local_c, local_d, area, half_distance) in specifications.items():
        condition = bc[face]
        if condition.kind in {"dirichlet", "inlet"}:
            rates[face] += float(np.sum(local_d * area / half_distance * (condition.value - local_c)))
        elif condition.kind == "robin":
            rates[face] += float(np.sum(condition.mass_transfer_m_s * area * (condition.value - local_c)))
    if velocity_x_m_s > 0:
        if bc["x_min"].kind == "inlet":
            rates["x_min"] += float(np.sum(velocity_x_m_s * grid.face_area_x_m2[:, :, 0] * bc["x_min"].value))
        if bc["x_max"].kind == "outlet":
            rates["x_max"] -= float(np.sum(velocity_x_m_s * grid.face_area_x_m2[:, :, -1] * c[:, :, -1]))
    return rates

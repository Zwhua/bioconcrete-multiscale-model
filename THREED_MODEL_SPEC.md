# Three-dimensional reactive-transport model specification

Status: frozen scientific specification for the v0.6.0 development cycle. The
model outputs are **uncalibrated 3D model output** and **not experimental data**.

## Scope and coordinate convention

The reference model resolves a crack in three spatial directions. `x` runs from
the exposed entrance toward the crack tip (or from inlet to outlet), `y` spans
the aperture between the two crack walls, and `z` runs along crack depth/front.
All public arrays use `(z, y, x)` and are flattened and restored in C order.
Files and run manifests must record `axis_order: [z, y, x]` and
`flatten_order: C`.

The model answers whether aperture-wise and depth-wise transport gradients
alter reaction, mineral distribution, closure, and internal connectivity. The
2.5D production model instead solves aperture-averaged fields on `(x,z)` while
retaining the local aperture `b(x,z,t)`. It is admissible only when aperture
mixing is fast relative to reaction, wetting and in-plane transport.

This is a reactive-transport and aperture-evolution model. It is not a complete
3D fracture-mechanics model, a single-crystal or single-cell agent model, a
general CT mesher, a full 3D Sobol engine, or a predictor of structural strength.

## Geometry and topology

The first implemented geometry is rectangular with dimensions taken only from
`transport.crack_length_mm`, `transport.crack_width_mm`, and
`transport.crack_depth_mm`. The legacy `out_of_plane_thickness_mm` has no role
in 3D volume. Cell volumes and face areas are physical SI quantities; inventories
are integrated as `sum(cell_value * cell_volume)`.

Two topologies have distinct semantics:

| Topology | x-min | x-max | Permitted hydraulic result |
| --- | --- | --- | --- |
| `blind_crack` | exposed entrance | no-flux tip | penetration depth, open volume, local cubic-law proxy only |
| `through_crack` | explicit inlet | explicit outlet | pressure-driven total flow and effective transmissivity |

A blind crack must never be assigned or reported a fictitious steady through-flow.

## Governing transport equation

For each mobile inventory (lactate, oxygen, calcium and total inorganic carbon
in equilibrium carbonate mode),

\[
\frac{\partial(\phi S_w C_i)}{\partial t}
+\nabla\cdot(\mathbf{u}C_i)
=\nabla\cdot(\phi S_w D_{i,\mathrm{eff}}\nabla C_i)+R_i.
\]

The finite-volume face diffusivity is the harmonic mean
`D_f = 2 D_P D_N / (D_P + D_N)` and the oriented diffusive face flux is
`F_PN = A_f D_f (C_N-C_P)/d_PN`. Physical face area, cell volume, centre
distance and boundary half-cell distance are mandatory. The initial advection
scheme is first-order upwind in `x`; diffusion uses a seven-point 3D stencil.

In equilibrium carbonate mode only total inorganic carbon is transported and
carbonate fractions are recomputed after transport. Kinetic species may be
transported independently only in explicit kinetic mode.

## Face boundary conditions

| Face | Default | Meaning |
| --- | --- | --- |
| `x_min` | exposed oxygen/DIC Robin boundary | exchange with the environment or explicit inlet |
| `x_max` | no flux for blind crack | sealed tip; outlet only for through crack |
| `y_min`, `y_max` | crack wall | no ordinary penetration flux; wall reaction/deposition may be accounted separately |
| `z_min`, `z_max` | no flux | finite crack-depth boundaries unless explicitly configured |

Supported conditions are `no_flux`, `dirichlet`, `robin`, `inlet`, `outlet`, and
`crack_wall`, configured per species. Oxygen Robin exchange obeys
`-D_O grad(C_O) dot n = k_L (C_O* - C_O)`. The 3D default is
`boundary_robin`; it must not be combined with legacy volumetric oxygen transfer.
`legacy_volumetric` is retained only for controlled dimensional regression.

## Reaction and splitting

All dimensions call one shared 22-state local reaction kernel. A 3D calculation
must not duplicate reaction equations or construct one global BDF system over
all voxels. Local cells are integrated in deterministic batches (default 64).
Production coupling uses Strang splitting: half transport, full local reaction,
half transport. Lie splitting remains available for legacy regression. Steps
end exactly at output, wet/dry transition, and checkpoint events.

Failed transport or reaction records time and cell indices, halves the step and
retries a bounded number of times. Persistent failure terminates with a failure
manifest; failed cells are never silently zeroed.

## Deposition, aperture and conservation

Calcite and C-S-H are separate inventories. In `surface_resolved` mode, for each
`(x,z)` column,

\[
V_s=\sum_y(V_{calcite}+V_{CSH}),\quad V_w=f_wV_s,\quad
V_{bulk}=V_s-V_w.
\]

Symmetric wall deposition gives `V_top = V_bottom = V_w/2` and
`b=max(b0-(V_top+V_bottom)/A_xz,b_min)`. The fixed-grid implementation retains
existing solute and solid mass after a column closes. It changes the local fluid
fraction, transport resistance and closed-column mask rather than deleting cells.

Required invariants are `V_wall + V_bulk = V_total`, independent calcite/C-S-H
accounting, carbon and calcium balance, and nonnegative finite states. Surface
closure is not equivalent to loss of internal connectivity. Overall closure uses
physical area weighting, never an unweighted array mean.

## Evidence boundary

Every 3D field, summary and figure carries model version, Git commit, config and
geometry hashes, grid resolution, and the labels `Uncalibrated 3D model output`
and `Not experimental data`. Synthetic rectangular or rough geometries are
numerical tests, not CT specimens. No structural-strength or experimental-
validation claim follows from this model alone.

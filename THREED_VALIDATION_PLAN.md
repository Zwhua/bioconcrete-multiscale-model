# 3D validation and verification plan

This plan is preregistered for v0.6.0. Passing numerical verification establishes
implementation credibility, not experimental validity.

## Ordered gates

1. **Gate A — specification:** coordinate, topology, face boundary conditions,
   2.5D/3D roles and exclusions are unambiguous.
2. **Gate B — reaction kernel:** all legacy tests and frozen 0D/1D/2D fixtures
   pass; single-cell and batched reaction agree; only one reaction source exists.
3. **Gate C — pure 3D transport:** manufactured solutions, boundary mass closure,
   real `z` diffusion and linear-solver convergence pass before reaction coupling.
4. **Gate D — coupled 3D:** dimensional reductions, conservation, nonnegativity,
   restart equivalence, grid convergence and time convergence pass before formal
   figures are produced.
5. **Gate E — 2.5D:** representative cases differ from 3D by less than 5% inside
   the declared applicability domain, warn outside it, and respect topology.
6. **Gate F — release candidate:** complete tests/CI, clean-commit result manifest,
   evidence-labelled figures, and no false validation language.

No later scientific layer is accepted when its preceding gate fails.

## Numerical tests and thresholds

| Test | Metric | Acceptance |
| --- | --- | --- |
| closed carbon balance | relative inventory residual | `< 0.5%` |
| closed calcium balance | relative inventory residual | `< 0.5%` |
| open-system balance | state change minus boundary and reaction integrals | `< 0.5%` or documented stricter scale-aware absolute tolerance |
| 3D diffusion MMS | L1, L2, Linfinity order | L2 observed order `>= 1.8` |
| upwind advection MMS | L1, L2, Linfinity order | observed order `>= 0.9` |
| pure transport reduction | 3D vs lower-dimensional field | `rtol <= 1e-6` |
| coupled reduction | summary relative error | `< 1%`, or convergence explanation |
| grid convergence | medium vs fine summary metrics | `< 5%` |
| time convergence | `dt/2` vs `dt/4` summary metrics | `< 5%` |
| 2.5D applicability | representative 2.5D vs 3D metrics | `< 5%`; `>10%` rejects automatic use |
| restart | continuous vs resumed fields | tolerance-level equality |

Grid refinement changes all directions. The nominal sequence displayed as
`nx x ny x nz` is `21x3x9`, `41x5x17`, `81x7x33`; internal shape is always
`(nz,ny,nx)`. Time refinement uses `dt`, `dt/2`, `dt/4`.

## Required test families

- constant and zero-gradient no-flux invariance;
- pure diffusion analytic solution and 3D manufactured solutions;
- nonuniform `z` data diffuse while uniform `z` data remain uniform;
- direct and iterative solver agreement with recorded convergence diagnostics;
- one voxel to 0D, `z`-uniform to 2D, and `y,z`-uniform to 1D;
- zero dose, spores, calcium, carbon, oxygen and wall-deposition limits;
- fully closed/zero-pressure hydraulic limits and topology-specific reporting;
- wall/bulk/total deposition partition and independent calcite/C-S-H closure;
- finite, nonnegative, no-NaN output; fixed-seed and serial/parallel repeatability;
- checkpoint hash rejection and continuous/resumed equivalence.

## Prospective experimental falsification

Micro-CT should quantify the depth-dependent mineral/void field and test whether
surface closure coexists with open interior paths. Controlled water-flow tests
on a through-crack specimen should measure pressure drop and flow before/after
healing. Dissolved oxygen or a validated oxygen-sensitive method should constrain
penetration. Capsule layouts (surface-concentrated, uniform-depth and layered)
must have the same independently measured total inventory. These measurements
can falsify predicted depth uniformity, connectivity and 2.5D applicability;
until performed, outputs remain uncalibrated predictions.

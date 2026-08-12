# V0.4.0 Scientific Audit

Evidence status: model and numerical evidence only. Team wet-lab observations: 0.
The default 28-day values below are uncalibrated model outputs.

## Findings

1. **The 2.0819% closure is almost entirely C-S-H filling.** The current 0D
   baseline assigns 2.08023 percentage points to C-S-H and 0.00168 percentage
   points to calcite. This is a model-prior decomposition, not measured phase
   attribution.
2. **Calcite is currently too small to affect total closure materially.** The
   predicted calcite mass is 0.03638 mg for the configured crack volume. Paired
   mineral mass and crack microscopy are required to test the precipitation and
   wall-deposition assumptions separately.
3. **Activity is not the dominant modeled limitation.** Increasing effective
   activity has little leverage because available inventory, release,
   precipitation, wall deposition and crack geometry constrain the downstream
   conversion to closure. The classifier reports the dominant limitation per
   scenario instead of forcing an activity-limited interpretation.
4. **The continuous-wetting conclusion is conditional.** Its modeled benefit
   depends on the oxygen boundary-transfer and wet-state diffusion assumptions.
   Dissolved oxygen and moisture time series are required to falsify this
   mechanism.
5. **More mineral mass with lower closure in a wider crack is geometrically
   consistent.** Closure is wall-deposit thickness divided by initial aperture;
   a larger void can contain more mass while closing by a smaller fraction.
6. **Fixed-total-inventory semantics are shared across scales.** Calcium
   lactate, C-S-H payload, spores and initial active units scale together, and
   initial integrated inventory is tested across 0D, 1D and 2D grids. Spatial
   predictions still require matched boundaries and discretization before
   quantitative cross-scale comparison.
7. **Transport mappings have theory but need material-specific validation.**
   Crack transmissivity uses the cubic law and matrix permeability uses a
   relative Kozeny-Carman relation. Neither mapping is a substitute for measured
   permeability or sorptivity in this material.
8. **Wall-deposition fraction is not currently identifiable from closure alone.**
   It is confounded with mineral volume and crack geometry. At minimum, matched
   CaCO3 mass and width maps at multiple times are needed.
9. **C-S-H and calcite are accounted separately.** Solid and wall volumes and
   their closure contributions are reported independently. An abiotic C-S-H
   control is still necessary to validate that the physical processes are not
   experimentally double-counted.
10. **The model has more uncertain parameters than independent observations.**
    Local Fisher diagnostics, parameter correlations and measurement
    recommendations are implemented, but scientific identifiability cannot be
    claimed until real observations with audited uncertainties are available.

## Design Consequence

The first useful experiment is not a broad activity screen. It is a compact
factorial comparison of complete repair material, no-C-S-H control and abiotic
C-S-H control at matched total inventory, with paired substrate, oxygen,
calcium, CaCO3 mass and crack-width observations. This directly tests the
current dominant closure attribution and reveals whether release, mineralization
or wall deposition should be improved next.

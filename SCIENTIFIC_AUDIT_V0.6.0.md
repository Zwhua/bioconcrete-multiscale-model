# Scientific audit — v0.6.0 development candidate

The implementation resolves `(x,y,z)` on cell-centred finite volumes with public
array order `(z,y,x)`. Geometry integrates physical volume, diffusion uses face
areas and harmonic face coefficients, and the test suite demonstrates a real
depth gradient rather than a copied 2D field. Calcite and C-S-H wall/bulk volumes
are separately retained. Blind cracks cannot report through-flow.

The legacy v0.5.1 0D/1D/2D artifacts are frozen at commit `1890526`; regression
tests protect their reported summaries. Local reactions use the existing single
scientific implementation through a dimension-independent batched entry point.

Gate D is passed by the versioned validation report: carbon/calcium balance,
dual dimensional-reduction adapters, nonnegativity, restart equivalence, and
grid/time convergence meet their frozen thresholds. Formal fields use
metadata-complete consolidated Zarr and evidence-labelled static figures.

Release remains blocked by Gate E/F work: the complete 2.5D reactive-transport
production solver, its accepted 2.5D--3D applicability table, clean-commit CI,
and release-candidate manifest review remain incomplete. No v0.6.0 version bump
is scientifically permitted yet.

All 3D development output is `uncalibrated 3D model output` and
`not experimental data`. It does not establish strength recovery, complete
hydraulic sealing, CT validation, or wet-lab validation.

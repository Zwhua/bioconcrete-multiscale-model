# v0.5.0 Release Candidate Checklist

## Implemented

- [x] Central evidence-state vocabulary and claim gating.
- [x] Counterfactual model-response bottleneck analysis.
- [x] Numerical Jacobian and greedy D-optimal experiment design.
- [x] Time- and condition-resolved six-structure comparison.
- [x] Anonymous design-category to parameter mapping.
- [x] Resumable multiprocessing sensitivity and design tasks.
- [x] V5 release manifest, static figures and evidence dashboard.
- [x] V0.2 figures removed from the current-results narrative.

## Required Before Scientific Result Release

- [ ] Verify and curate a licensed public calibration workbook.
- [ ] Complete specimen-grouped public calibration and freeze parameters.
- [ ] Complete independent external evaluation without refitting.
- [ ] Run 1024-base-sample Sobol analysis to completion.
- [ ] Run all 1,728 preregistered scenarios to completion.
- [ ] Complete prior predictive, counterfactual and numerical D-optimal runs.
- [ ] Regenerate Figures 2-8 from completed V5 artifacts.
- [ ] Run full regression, physical validation and CI on a clean commit.

Unchecked items block calibrated or externally evaluated claims. They do not
permit smoke outputs or synthetic fixtures to be promoted as scientific data.

## Git Release

The software release may be committed and tagged while the unchecked scientific
analyses remain explicitly incomplete. Their smoke outputs must not be promoted
as formal evidence.

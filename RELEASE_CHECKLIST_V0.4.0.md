# v0.4.0 Release Checklist

## Completed

- [x] Fixed-total-inventory dose semantics across 0D, 1D and 2D.
- [x] Calcite and C-S-H closure contributions separated.
- [x] Physical limit and deterministic-seed tests.
- [x] Local sensitivity, Fisher information and correlation diagnostics.
- [x] Estimability labels and measurement recommendations.
- [x] Resumable prior predictive uncertainty with epistemic/scenario separation.
- [x] Six model structures with data-gated AIC/AICc.
- [x] Seven-class bottleneck analysis.
- [x] D-optimal experiment ranking and decision-support tables.
- [x] Deterministic, resumable and multiprocessing design matrix.
- [x] Machine-readable run manifests and dirty-worktree warning.
- [x] Python 3.8/3.11/3.12 CI plus manual formal-analysis workflow.
- [x] Prospective predictions and prospective DBTL wording.
- [x] 62 local tests, compile check, whitespace check and quick validation pass.
- [x] Package, citation, README and documentation set to `0.4.0`.

## Evidence Still Required

- [ ] Curate and license-check public calibration dataset A.
- [ ] Freeze fitted shared parameters after specimen-grouped internal testing.
- [ ] Run independent dataset B without refitting.
- [ ] Fit the image measurement-error model with krkCMd.
- [ ] Cross-check the analytical chemistry surrogate against PHREEQC.
- [ ] Execute 1024-base-sample formal Sobol convergence analysis.
- [ ] Execute all 1,728 preregistered scenarios on a frozen calibrated config.
- [ ] Execute future team wet-lab measurements; current team observations: 0.

Unchecked items are not release blockers for the software interfaces, but they
are blockers for calibration, external-validation and experimental-performance
claims.

## Release Commands

```powershell
python -m compileall -q bioconcrete tests
python -m unittest discover -s tests -v
python -m bioconcrete validate --config model_config.json
git diff --check
git status --short
git tag -a v0.4.0 -m "Release v0.4.0"
git push origin main
git push origin v0.4.0
```

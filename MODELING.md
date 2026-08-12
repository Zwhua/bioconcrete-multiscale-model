# BioConcrete multiscale model

This package implements a physicochemical self-healing concrete model. Its
biological boundary is an anonymous population-scale activity interface; no
biological design information or construction details are stored or inferred.

## Model layers

- `0d`: capsule release, activity, aerobic calcium-lactate utilization,
  carbonate speciation, Portlandite dissolution, calcite precipitation, and
  material repair metrics.
- `1d`: implicit finite-volume transport from the crack mouth to its tip,
  coupled to local BDF reaction solves.
- `2d`: a true length-by-width finite-volume grid with discrete capsule source
  fields. It is not an interpolation of the 1D result.

The tracked calcium-lactate pathway is ammonia-free. `ammonium_mol_m3` is kept
as a diagnostic state and must remain exactly zero.

## Quick start

```powershell
python -m bioconcrete default-config --output model_config.json
python -m bioconcrete prepare-data
python -m bioconcrete build-geochem-grid
python -m bioconcrete simulate --level 0d --config model_config.json
python -m bioconcrete simulate --level 1d --config model_config.json
python -m bioconcrete simulate --level 2d --config model_config.json
python -m bioconcrete validate --config model_config.json
```

Simulation folders contain `state.csv`, `summary.json`, `diagnostics.json`, the
exact configuration, a figure, and `REPORT.md`.

## Public calibration and external validation

The legacy CSV interface remains available. The v0.2 evidence workflow uses the
traceable public observation schema documented in `data/public/README.md` and
keeps project experiments, public calibration, and external validation separate.

```text
lactate_mM, cfu_mL, caco3_mg, healing_ratio,
permeability_ratio, sorptivity
```

Run calibration with:

```powershell
python -m bioconcrete calibrate --data experiment.csv --bootstrap 20
```

Public database values are parameter priors. The software does not treat them
as measurements of the project material. An 80% healing ratio is shown only as
an evaluation target and is never hard-coded into the equations.

`validate` checks mass balance, limiting cases, time-step convergence, and grid
convergence. Add `--full` for the expensive full-duration refined-grid check.

## Geochemical backend

`build-phreeqc-grid` records the backend that actually generated the table. On
the current machine no PHREEQC executable is installed, so the table is labelled
`analytical_surrogate`. `compare-geochem-backends` therefore reports
`claim_allowed: false`; no PHREEQC-coupling claim is made.

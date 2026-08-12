# Model

## Evidence boundary

This repository has no team wet-lab time series. It separates five evidence
classes: project experiment (currently empty), public calibration data,
independent external validation data, literature priors, and model predictions.
Synthetic records are used only in software tests and are never reported as
scientific calibration evidence.

## What the model predicts

The model links capsule release, anonymous effective activity, carbonate
chemistry, calcite precipitation, transport, and crack-wall deposition. It
provides 0D, 1D, and true 2D solvers. Its principal outputs are solid fill,
wall deposition thickness, crack closure, permeability, transmissivity, and
total mineral mass.

The activity interface contains only population-scale aggregate parameters:
effective turnover, effective half-saturation, active-unit concentration,
activity multiplier, response delay, basal leakage, and activation duration.
It represents material-scale effective response only and cannot identify an
underlying implementation.

## Dynamic environment

Water activity, oxygen and pH drive an environment signal. The model tracks
signal duration, response delay, activation state, memory and cumulative
resource consumption. The default `AND` mode requires water suitability,
oxygen-rise and pH-drop gates simultaneously; `OR` and `static_suitability` are
retained only as explicit ablation modes. `true_activation_index` and
`false_activation_index` are both evaluated against the same strict AND
reference. The module is a phenomenological activation surrogate, not a
validated model of a particular biological circuit. Dynamic pH is solved from dissolved inorganic carbon and
total alkalinity after every accepted reaction step. Fixed pH remains an
explicit comparison mode.

## Geometry

Mineral concentration is converted to total mass using crack length, aperture
and the appropriate unresolved depth or out-of-plane thickness. Repair no
longer equates solid volume fraction with closure. The wall-deposited portion
gives a thickness

```text
delta_one_wall = wall deposit volume / total opposing-wall area
closure = clip(2 * delta / initial aperture, 0, 1)
```

Here `total opposing-wall area` is the sum of both crack faces and `delta` is
the thickness on one face. The factor of two therefore appears only in the
closure equation. Wall and non-wall solid volumes are reported separately and
sum to the total deposited-solid volume.

`healing_ratio` remains only as a deprecated file-format alias.

## Cross-scale accounting

The configured dose is independent of the number of discrete capsule sources.
Source count changes spatial heterogeneity, while each profile is normalized so
its volume-weighted mean inventory is unchanged by grid refinement. Total
calcite, carbon, calcium and capsule inventories are computed as sums of local
concentration times cell volume. Comparisons across 0D, 1D and 2D require the
same physical crack volume; in particular, the 2D out-of-plane thickness must
match the depth used by 0D/1D.

The v0.4 model candidate also fixes total repair-agent inventory when crack
geometry changes. One dosage multiplier scales substrate, releasable C-S-H,
spores and initial active material together. Calcite and C-S-H closure
contributions are reported separately; C-S-H-only filling is explicitly
nonbiological model output.

## Calibration and blind validation

Dataset A is split by specimen, so time points from one specimen cannot cross
between training and internal testing. Shared parameters are frozen with a
SHA-256 digest. Dataset B can only enter the external-validation command, which
contains no optimizer and rejects modified frozen configurations. It reports
the mechanistic model beside zero-mineralization and first-order baselines.

Dataset C is restricted to crack-width measurement error. Its uncertainty is
not used to fit reaction kinetics.

## Uncertainty and design

Formal sensitivity uses SALib Morris trajectories and Saltelli/Sobol sampling,
with bootstrap confidence intervals and explicit non-convergence flags. The
preregistered matrix fixes 1,728 scenarios before external validation. Pareto
labels are prospective design suggestions, not verified experimental choices.

## Current limitation

Zenodo returned HTTP 403 in the automated environment, so no calibration or
external-validation score is currently claimed. The complete data pipeline and
manual download paths are provided in `data/public/MANUAL_DOWNLOADS.md`.
The geochemical table currently uses `analytical_surrogate`; without an actual
PHREEQC executable, reports must not describe it as PHREEQC coupled.

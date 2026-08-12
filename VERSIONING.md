# Versioning Policy

This project follows Semantic Versioning using `MAJOR.MINOR.PATCH`.

- `PATCH`: backward-compatible fixes, documentation corrections, tests, and
  small maintenance changes that do not alter scientific interpretation.
- `MINOR`: backward-compatible model capabilities, new outputs or commands,
  calibration workflow changes, and scientific behavior changes that require
  users to review regenerated results.
- `MAJOR`: incompatible configuration, CLI, data-schema or model-contract
  changes that require user migration.

Every released change must update the same version in:

- `pyproject.toml`;
- `bioconcrete/__init__.py`;
- `CITATION.cff`;
- the README version badge and current-version statement.

Release tags use the `vMAJOR.MINOR.PATCH` form. Development commits may be made
between releases, but the version is incremented before those changes are
pushed as a declared release. Historical result labels retain the version that
actually produced them.

## Current Release

`v0.4.0` adds fixed-total-inventory dose semantics, practical identifiability,
prior predictive uncertainty, mechanistic model comparison, bottleneck
classification, D-optimal experiment ranking, resumable parallel design
matrices, unified run manifests, and prospective decision-support artifacts.
It does not claim that public calibration, external validation, or a team wet-lab
DBTL cycle has completed.

`v0.3.0` remains the historical release that introduced multi-output staged
calibration, explicit environmental gate logic, volume-conserving deposition
geometry, and reviewed-public-data requirements.

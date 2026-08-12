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

`v0.3.0` formalizes the multi-output staged calibration workflow, explicit
environmental gate logic, volume-conserving deposition geometry, cross-scale
inventory accounting, and reviewed-public-data requirements. It does not claim
that public calibration or external validation has completed; those stages
remain blocked until the selected source files are downloaded and curated.

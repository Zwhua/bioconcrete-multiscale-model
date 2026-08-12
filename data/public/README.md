# Public evidence data

`DATASETS.yml` separates the evidence roles before any model fitting:

- `transet_18clsu02`: public calibration data;
- `marine_external`: external validation data, never used for fitting shared parameters;
- `krkcmd`: measurement-error data only.

Raw downloads, extraction workspaces, receipts, and derived observations are
ignored by Git. Run `fetch-public-data` to create SHA-256 receipts and
`prepare-public-data` to create traceable local tables. Every normalized row
keeps its source file, sheet, and row number. No file is described as a project
experiment: this project currently has no team wet-lab observations.

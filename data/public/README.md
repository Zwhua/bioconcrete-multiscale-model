# Public evidence data

`DATASETS.yml` separates the evidence roles before any model fitting:

- `bath_aea_spores`: preferred public calibration candidate; repository metadata,
  license and workbook fields must be verified before calibration;
- `transet_18clsu02`: public calibration data;
- `marine_external`: external validation data, never used for fitting shared parameters;
- `krkcmd`: measurement-error data only.

Raw downloads, extraction workspaces, receipts, and derived observations are
ignored by Git. Run `fetch-public-data` to create SHA-256 receipts and
`prepare-public-data` to create traceable candidate tables. Generic extraction
is discovery-only and is marked `candidate_only`; `calibrate-public` rejects it.
Formal use requires manual review under [DATA_CURATING_PROTOCOL.md](../../DATA_CURATING_PROTOCOL.md)
and an approved curated data dictionary. Every normalized row
keeps its source file, sheet, and row number. No file is described as a project
experiment: this project currently has no team wet-lab observations.

The normalized observation table supports `lactate_mM`, `calcite_mass_mg`,
`crack_closure_ratio` (or initial/current crack widths), `permeability_ratio`,
`ph`, `activation_state`, and `cumulative_activity_h`. Missing values remain
empty; the preparation pipeline must not invent measurements. Per-output
standard deviations use `<output>_sd`. The generic `measurement_sd` field is a
crack-width uncertainty in mm and is propagated to closure uncertainty using
the initial width.

Current acquisition status: no selected public raw file is present locally, so
public calibration and external validation are blocked. Empty directories are
not evidence, and no synthetic replacement dataset is generated.

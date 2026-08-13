# Manual public-data download

The current automated environment cannot reliably complete repository HTTPS
downloads. This is an access limitation, not evidence that the records are
absent. Download in a normal browser and place the files at the paths below:

1. Preferred calibration candidate: open
   <https://researchdata.bath.ac.uk/1087/>. Verify the displayed DOI, license,
   workbook filename and checksum before placing the spreadsheet under
   `data/public/raw/bath_aea_spores/`. Required fields are specimen identity,
   initial crack width, healing time, crack healing ratio, treatment/control
   group and water penetration recovery. Do not approve it for calibration
   until these fields and redistribution terms have been checked.
2. Fallback calibration candidate: open <https://zenodo.org/records/3471960> and download
   `TranSET_18CLSU02_Data.zip` to
   `data/public/raw/transet_18clsu02/TranSET_18CLSU02_Data.zip`.
3. Dataset B: open <https://zenodo.org/records/11305154> and download
   `02_Mar_Microscope measurements.xlsx` and `06_Mar_Mass change.xlsx` to
   `data/public/raw/marine_external/`.
4. Measurement dataset: open <https://zenodo.org/records/14568863>. Access may
   be restricted. Download only the tabular crack-width labels, when permitted,
   to `data/public/raw/krkcmd/`. Do not download images for this model task.

Then run candidate discovery:

```powershell
python -m bioconcrete prepare-public-data --dataset transet_18clsu02
python -m bioconcrete prepare-public-data --dataset marine_external
```

Before scientific use, compare local SHA-256 hashes with a generated receipt or
the publisher's checksums. Raw files and local derived tables remain ignored by
Git. The generated rows are marked `candidate_only` and cannot be calibrated.
Review every `source_location`, complete the curated observations and dictionary
under `data/public/curated/`, and set dictionary `reviewer_status=approved` only
after independent source verification. The conservative generic extractor
intentionally leaves ambiguous fields empty.

Current audit (2026-08-12): the raw folders contain no selected files. Automated
metadata access was intermittent; the selected file endpoint returned HTTP 403
and the record page later returned HTTP 429. Calibration freeze and external
validation therefore remain blocked, and no replacement observations were made.

# Manual public-data download

The current automated environment receives HTTP 403 from Zenodo. This is an
access limitation, not evidence that the records are absent. Download in a
normal browser and place the files at the exact paths below:

1. Dataset A: open <https://zenodo.org/records/3471960> and download
   `TranSET_18CLSU02_Data.zip` to
   `data/public/raw/transet_18clsu02/TranSET_18CLSU02_Data.zip`.
2. Dataset B: open <https://zenodo.org/records/11305154> and download
   `02_Mar_Microscope measurements.xlsx` and `06_Mar_Mass change.xlsx` to
   `data/public/raw/marine_external/`.
3. Measurement dataset: open <https://zenodo.org/records/14568863>. Access may
   be restricted. Download only the tabular crack-width labels, when permitted,
   to `data/public/raw/krkcmd/`. Do not download images for this model task.

Then run:

```powershell
python -m bioconcrete prepare-public-data --dataset transet_18clsu02
python -m bioconcrete prepare-public-data --dataset marine_external
```

Before scientific use, compare local SHA-256 hashes with a generated receipt or
the publisher's checksums. Raw files and local derived tables remain ignored by
Git. Review every `source_location` row before calibration; the conservative
generic extractor intentionally leaves ambiguous fields empty.

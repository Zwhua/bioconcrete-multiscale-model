# Public Data Curation Protocol

Status: preregistered before inspection of external-validation outcomes.

## Evidence Roles

- `transet_18clsu02`: public calibration and specimen-grouped internal testing.
- `marine_external`: external validation only; shared parameters remain frozen.
- `krkcmd`: crack-width measurement-error estimation only.
- Project wet-lab time-series data: none at the time of this protocol.

Generic spreadsheet extraction produces candidate rows only. It is not a
scientific dataset and cannot enter calibration until a reviewer approves each
included field in the curated dictionary.

## Inclusion Criteria

Include a measurement only when the source identifies the specimen or an
unambiguous replicate, time point, endpoint definition, unit and treatment
condition. Initial and current crack widths must refer to the same measurement
location or to a clearly documented aggregate from the same specimen.

## Exclusion Criteria

Exclude values with ambiguous units, pooled groups without sample size,
unrecoverable figure axes, duplicate publication of the same observation,
missing specimen linkage, or conditions that cannot be distinguished from a
control. Record every exclusion and reason in the dictionary.

## Definitions

`initial_crack_width_mm` is the width at the start of the documented healing
period. `current_crack_width_mm` is measured at `time_d` using the source's
stated method. Closure is

```text
crack_closure_ratio = 1 - current_crack_width_mm / initial_crack_width_mm
```

Widths must be positive. Closure outside [0, 1] is retained in the raw review
sheet, flagged, and excluded from calibration unless the source documents a
valid interpretation.

## Replicates And Missingness

Keep biological/material replicates as separate specimens. Repeated locations
within one specimen may be retained as separate measurement rows only when the
location identifier is traceable; otherwise use the published aggregate and
its reported uncertainty. Do not impute missing scientific observations.

## Outliers

Do not remove observations solely because they fit poorly. Exclude only for a
predefined data-quality reason, such as transcription error, impossible unit,
or a source-marked failed specimen. Preserve the original value and reason.

## Units And Digitization

Normalize time to days, widths to mm, mass to mg, concentration to mM, and
ratios to fractions. Every conversion requires an explicit formula in the data
dictionary. Figure digitization must record figure/panel, software, axis scale,
pixel or point uncertainty, and an independent reviewer check.

## Splitting And Fitting

Calibration/internal-test assignment is grouped by `specimen_id`; no time point
from one specimen may cross splits. Stage-specific observations determine which
parameters may be fitted. Missing calcite mass fixes precipitation rate to its
prior. External-validation rows are never used for fitting or model selection.

## Freeze Record

Before external validation, record this file's SHA-256, curated files and
checksums, dictionary reviewer status, split mapping, software commit and frozen
configuration hash. Any later change creates a new protocol version.

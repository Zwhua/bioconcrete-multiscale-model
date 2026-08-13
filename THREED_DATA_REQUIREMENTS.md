# 3D data requirements and evidence policy

## Minimum calibration data

| Quantity | Minimum metadata | Intended constraint |
| --- | --- | --- |
| crack geometry / micro-CT | voxel size, orientation, segmentation method, specimen ID, time | aperture and internal connectivity |
| repair-agent placement | capsule coordinates/distribution, total mass or moles, batch | fixed total inventory and spatial profile |
| fluid boundary history | inlet concentration, pressure/flow, temperature, wet/dry timing | transport boundary conditions |
| oxygen and inorganic carbon | method, position/depth, time, uncertainty | Robin exchange and penetration |
| calcium/lactate | method, time, location, uncertainty | source/reaction inventories |
| mineral phase | calcite and C-S-H measured separately where possible | deposition partition |
| surface closure | registered image scale, threshold procedure, uncertainty | surface metric only |
| hydraulic response | inlet/outlet definition, pressure difference, flow, viscosity | through-crack transmissivity |

Rows need units, specimen/group identifiers, replicate identifiers, timestamps,
measurement uncertainty, preprocessing record, and provenance/licence. Training,
validation and test specimens must be separated by specimen rather than by row.

## Geometry inputs

`rectangular` is a model geometry. `aperture_field` accepts a validated 2D field
with physical units, axis convention and missing-value policy. `voxel_ct` is only
an input/validation interface until actual CT data meeting the metadata above are
provided. Analytic roughness is explicitly labelled synthetic and may only serve
verification or sensitivity analysis.

## Evidence classes

- `model_prior`: literature/database-derived parameter range;
- `project_hypothesis`: uncalibrated project choice;
- `prospective uncalibrated 3D model prediction`: frozen falsifiable prediction;
- `team_wet_lab`: admissible only for traceable measured rows;
- `public_calibration`: admissible only after frozen, reproducible calibration.

Missing data are reported, never fabricated or imputed as experimental evidence.
All current 3D outputs must state `uncalibrated 3D model output` and
`not experimental data`. Surface photographs alone cannot validate interior
closure, connectivity, phase identity or flow resistance.

## Recommended experiment

Use matched through-crack specimens with equal total capsule inventory and three
depth layouts: surface-concentrated, depth-uniform and layered. Acquire baseline
and post-healing micro-CT, registered surface images, and pressure-controlled flow;
record wet/dry and chemical boundary histories. The paired observations test the
specific 3D prediction that high surface closure need not imply low interior
connectivity and determine whether layout changes minimum depth closure.

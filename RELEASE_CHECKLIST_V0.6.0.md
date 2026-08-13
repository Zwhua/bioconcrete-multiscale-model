# v0.6.0 release checklist

Status: **not release-ready**. Development commits may be pushed for review, but
no v0.6.0 tag or GitHub Release is authorized until every release gate passes.

- [x] Gate A: 3D specification frozen
- [x] Gate B: shared reaction schema/kernel entry and v0.5.1 regression
- [x] Gate C: conservative 3D transport, real z gradient, flux closure, diffusion MMS
- [x] Gate D: complete dimensional reduction, conservation, restart and grid/time convergence
- [ ] Gate E: complete 2.5D reactive-transport solver and validated applicability warning
- [x] Required Xarray/Zarr field store and metadata
- [ ] Full representative validation from a clean commit
- [x] Required static evidence-labelled figures and optional interactive fallback
- [ ] Python 3.8/3.11/3.12 CI matrix and manual 3D analysis workflow
- [ ] Public/team wet-lab calibration (absence must remain explicit)
- [ ] Version synchronized to 0.6.0 only after every release gate passes

Current evidence is numerical verification of an uncalibrated development model,
not experimental validation and not a structural-mechanics prediction.

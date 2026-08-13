"""Evidence-labelled static and optional interactive 3D visualizations."""
from __future__ import annotations
import json
from pathlib import Path
from typing import Dict
import numpy as np

from .io_3d import require_gate_d


def _load(run: Path):
    try:
        import xarray as xr
    except ImportError as error:
        raise RuntimeError("render-3d requires the 'three-d' extra") from error
    run=Path(run)
    return xr.open_zarr(str(run/"fields.zarr"),consolidated=True)


def _label(dataset) -> str:
    return ("Model v0.6.0-development | Uncalibrated 3D model output | Not experimental data\n"
            "Git/config provenance in run_manifest.json | geometry {} | grid {}"
            .format(dataset.attrs.get("geometry_hash","unknown")[:12],
                    dataset.attrs.get("grid_shape","unknown")))


def render_static_3d(run: Path, output: Path = None, aperture_exaggeration: float=100.0) -> Dict[str,str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    run=Path(run); output=Path(output or run/"figures"); output.mkdir(parents=True,exist_ok=True)
    ds=_load(run); products={}
    fields=("oxygen_mol_m3","lactate_mol_m3","calcium_mol_m3",
            "calcite_mol_m3","csh_volume_fraction")
    for name in fields:
        a=np.asarray(ds[name].isel(time=-1));k,j,i=(v//2 for v in a.shape)
        fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
        panels=((a[k],"x-y"),(a[:,j,:],"x-z"),(a[:,:,i],"y-z"))
        vmin,vmax=float(np.nanmin(a)),float(np.nanmax(a))
        if vmax<=vmin:vmax=vmin+1e-30
        for ax,(panel,title) in zip(axes,panels):
            image=ax.imshow(panel,origin="lower",aspect="auto",vmin=vmin,vmax=vmax)
            ax.set_title(title);fig.colorbar(image,ax=ax,shrink=.75)
        fig.suptitle(name+"\n"+_label(ds),fontsize=8)
        for extension in ("png","svg"):
            path=output/("slices_{}.{}".format(name,extension));fig.savefig(str(path),dpi=180);products[path.name]=str(path)
        plt.close(fig)
    aperture=np.asarray(ds["aperture_m"].isel(time=-1));closure=np.asarray(ds["closure_ratio"].isel(time=-1))
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,data,title in ((axes[0],aperture*1e3,"Local aperture (mm)"),(axes[1],closure,"Closure ratio")):
        image=ax.imshow(data,origin="lower",aspect="auto");fig.colorbar(image,ax=ax);ax.set_title(title)
    fig.suptitle(_label(ds)+"\nAperture direction exaggerated {}x for visualization".format(aperture_exaggeration),fontsize=8)
    for extension in ("png","svg"):
        path=output/("aperture_closure."+extension);fig.savefig(str(path),dpi=180);products[path.name]=str(path)
    plt.close(fig)
    values=np.asarray(ds["closure_ratio"]);vmin=float(np.nanmin(values));vmax=float(np.nanmax(values))
    if vmax<=vmin:vmax=vmin+1e-30
    count=values.shape[0];fig,axes=plt.subplots(1,count,figsize=(4*count,4),squeeze=False,constrained_layout=True)
    for index,ax in enumerate(axes[0]):
        image=ax.imshow(values[index],origin="lower",aspect="auto",vmin=vmin,vmax=vmax)
        ax.set_title("t={:.3g} d".format(float(ds.time[index])));fig.colorbar(image,ax=ax,shrink=.7)
    fig.suptitle("Closure time evolution (fixed color scale)\n"+_label(ds),fontsize=8)
    for extension in ("png","svg"):
        path=output/("time_evolution."+extension);fig.savefig(str(path),dpi=180);products[path.name]=str(path)
    plt.close(fig)
    surface=closure[:,0] if closure.ndim==2 else closure
    open_path=np.asarray(ds["sealed_mask"].isel(time=-1))==0
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    axes[0].plot(np.ravel(surface));axes[0].set_title("Entrance/surface closure")
    axes[1].imshow(open_path,origin="lower",aspect="auto",cmap="binary");axes[1].set_title("Interior open columns")
    fig.suptitle(_label(ds),fontsize=8)
    for extension in ("png","svg"):
        path=output/("surface_vs_internal."+extension);fig.savefig(str(path),dpi=180);products[path.name]=str(path)
    plt.close(fig);ds.close();return products


def render_reduction_comparison(validation_report: Path, output: Path) -> Dict[str,str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    report=json.loads(Path(validation_report).read_text(encoding="utf-8"))
    comparison=report.get("grid_comparison_medium_fine",{})
    names=list(comparison);errors=[comparison[name].get("relative_error") or 0 for name in names]
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True);ax.bar(np.arange(len(names)),errors)
    ax.axhline(.05,color="red",linestyle="--",label="5% acceptance");ax.set_xticks(np.arange(len(names)))
    ax.set_xticklabels(names,rotation=25,ha="right");ax.set_ylabel("Relative error");ax.legend()
    ax.set_title("3D grid convergence / reduction evidence\nModel v0.6.0-development | Uncalibrated 3D model output | Not experimental data")
    products={}
    for extension in ("png","svg"):
        path=Path(output)/("reduction_comparison."+extension);fig.savefig(str(path),dpi=180);products[path.name]=str(path)
    plt.close(fig);return products


def render_pyvista_3d(run: Path, output: Path = None) -> Dict[str,str]:
    try:
        import pyvista as pv
    except ImportError:
        return {"pyvista":"not installed; static Matplotlib products are complete"}
    run=Path(run);output=Path(output or run/"figures");output.mkdir(parents=True,exist_ok=True)
    ds=_load(run);values=np.asarray(ds["calcite_mol_m3"].isel(time=-1))
    if not np.isfinite(values).any() or float(np.nanmax(values)-np.nanmin(values)) <= 1e-30:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        png=output/"calcite_isosurface.png";html=output/"calcite_isosurface.html"
        fig,ax=plt.subplots(figsize=(8,4));ax.axis("off")
        ax.text(.5,.55,"No non-trivial CaCO3 isosurface at this output time",ha="center")
        ax.text(.5,.35,"Model v0.6.0-development\nUncalibrated 3D model output\nNot experimental data",ha="center")
        fig.savefig(str(png),dpi=180);plt.close(fig)
        html.write_text("<html><body><h2>No non-trivial CaCO3 isosurface</h2>"
                        "<p>Model v0.6.0-development; Uncalibrated 3D model output; "
                        "Not experimental data.</p></body></html>",encoding="utf-8")
        ds.close();return {png.name:str(png),html.name:str(html),
                           "pyvista":"skipped empty isosurface"}
    grid=pv.ImageData(dimensions=np.array(values.shape[::-1])+1)
    grid.cell_data["calcite_mol_m3"]=values.ravel(order="C")
    threshold=float(np.nanpercentile(values,75));surface=grid.cell_data_to_point_data().contour([threshold],"calcite_mol_m3")
    if surface.n_points == 0:
        ds.close();return {"pyvista":"empty contour; static slice products remain authoritative"}
    plotter=pv.Plotter(off_screen=True);plotter.add_mesh(surface,scalars="calcite_mol_m3")
    plotter.add_text("Model v0.6.0-development\nUncalibrated 3D model output\nNot experimental data",font_size=9)
    png=output/"calcite_isosurface.png";html=output/"calcite_isosurface.html"
    plotter.show(screenshot=str(png),auto_close=False)
    products={png.name:str(png)}
    try: plotter.export_html(str(html));products[html.name]=str(html)
    except Exception as error: products["interactive_html"]="unavailable: {}".format(error)
    plotter.close();ds.close();return products


def render_formal_3d(run: Path, validation_report: Path, output: Path=None) -> Dict[str,str]:
    require_gate_d(validation_report)
    target=Path(output or Path(run)/"figures")
    products=render_static_3d(run,target);products.update(render_reduction_comparison(validation_report,target))
    products.update(render_pyvista_3d(run,target));return products

"""Atomic checkpoints and non-CSV field storage for 3D runs."""
from __future__ import annotations
import json, os
from pathlib import Path
from datetime import datetime
import shutil
import subprocess
import uuid
from typing import Any, Dict, Tuple
import numpy as np
from .config import ModelConfig
from .grid_3d import StructuredGrid3D
from .state import S, STATE_NAMES
from .chemistry import ph_from_alkalinity

def config_hash(config: ModelConfig) -> str:
    import hashlib
    payload=json.dumps(config.to_dict(),sort_keys=True,separators=(',',':')).encode()
    return hashlib.sha256(payload).hexdigest()

def write_checkpoint(path: Path, time_s: float, state: np.ndarray,
                     config: ModelConfig, grid: StructuredGrid3D,
                     ledger: Dict[str, Any] = None, step: int = 0,
                     output_index: int = 0, retry_count: int = 0) -> Path:
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    metadata=json.dumps({'time_s':time_s,'config_hash':config_hash(config),
                         'geometry_hash':grid.geometry_hash,'ledger':ledger or {},
                         'step':step,'output_index':output_index,
                         'retry_count':retry_count})
    with temporary.open('wb') as handle:
        np.savez_compressed(handle,state=state,metadata=np.array(metadata))
        handle.flush();os.fsync(handle.fileno())
    temporary.replace(path);return path

def read_checkpoint(path: Path, config: ModelConfig, grid: StructuredGrid3D,
                    full: bool = False):
    with np.load(Path(path),allow_pickle=False) as data:
        metadata=json.loads(str(data['metadata']));state=np.array(data['state'])
    if metadata['config_hash'] != config_hash(config): raise ValueError('checkpoint config hash mismatch')
    if metadata['geometry_hash'] != grid.geometry_hash: raise ValueError('checkpoint geometry hash mismatch')
    if state.shape != (grid.size,state.shape[1]): raise ValueError('checkpoint state shape mismatch')
    if full:
        return float(metadata['time_s']),state,metadata
    return float(metadata['time_s']),state

def estimate_peak_memory_gb(grid: StructuredGrid3D, saved_times: int,
                            state_count: int=22, save_full_state: bool=True) -> float:
    cells=grid.size; state=cells*state_count*8
    sparse_matrix=cells*7*(8+4)+cells*4
    saved=state*saved_times if save_full_state else 0
    work=state*4+sparse_matrix*2
    return float((state+saved+work)/1024**3)


FIELD_METADATA = {
    "oxygen_mol_m3": ("mol m-3", "Dissolved oxygen", "mobile"),
    "lactate_mol_m3": ("mol m-3", "Dissolved lactate", "mobile"),
    "calcium_mol_m3": ("mol m-3", "Dissolved calcium", "mobile"),
    "inorganic_carbon_mol_m3": ("mol m-3", "Total inorganic carbon", "mobile"),
    "calcite_mol_m3": ("mol m-3", "Calcite concentration", "solid"),
    "csh_volume_fraction": ("m3 m-3", "C-S-H volume fraction", "solid"),
}


def _optional_storage_modules():
    try:
        import xarray as xr
        from numcodecs import Blosc
        import zarr  # noqa: F401
    except ImportError as error:
        raise RuntimeError(
            "Formal 3D storage requires the 'three-d' extra: "
            "pip install '.[three-d]'. CSV/NPZ fallback is intentionally disabled."
        ) from error
    return xr, Blosc


def result_to_dataset(result):
    """Convert a simulation result to a metadata-complete Xarray Dataset."""
    xr, _ = _optional_storage_modules()
    grid = result.geometry
    coordinates = {"time": result.times_d, "z": grid.z_m, "y": grid.y_m, "x": grid.x_m}
    variables = {}
    for name, (unit, long_name, state_class) in FIELD_METADATA.items():
        variables[name] = (("time", "z", "y", "x"), result.state[..., S[name]], {
            "unit": unit, "long_name": long_name, "state_class": state_class,
            "evidence_class": "uncalibrated_model_output",
        })
    carbon=result.state[...,S["inorganic_carbon_mol_m3"]]
    alkalinity=result.state[...,S["total_alkalinity_mol_m3"]]
    ph=ph_from_alkalinity(carbon.ravel(),alkalinity.ravel(),
        result.config.environment.temperature_c,result.config.environment.ph_minimum,
        result.config.environment.ph_maximum,strict=False).reshape(carbon.shape)
    variables["ph"] = (("time","z","y","x"),ph,{"unit":"1","long_name":"Charge-balanced pH",
        "state_class":"derived_chemistry","evidence_class":"uncalibrated_model_output"})
    calcite_fraction = (result.state[..., S["calcite_mol_m3"]]
        * result.config.chemistry.calcite_molar_mass_kg_mol
        / result.config.chemistry.calcite_density_kg_m3)
    solid = np.clip(calcite_fraction + result.state[..., S["csh_volume_fraction"]], 0, 1)
    porosity = np.maximum(result.config.transport.porosity_initial * (1-solid),
                          result.config.transport.porosity_minimum)
    variables["porosity"] = (("time","z","y","x"), porosity, {
        "unit":"1", "long_name":"Effective porosity", "state_class":"geometry_proxy",
        "evidence_class":"uncalibrated_model_output"})
    variables["aperture_m"] = (("time","z","x"), result.aperture_history_m, {
        "unit":"m", "long_name":"Local crack aperture", "state_class":"geometry",
        "evidence_class":"uncalibrated_model_output"})
    variables["closure_ratio"] = (("time","z","x"), result.closure_history, {
        "unit":"1", "long_name":"Local aperture closure", "state_class":"geometry",
        "evidence_class":"uncalibrated_model_output"})
    variables["sealed_mask"] = (("time","z","x"), result.aperture_history_m <= 1e-9, {
        "unit":"1", "long_name":"Closed-column mask", "state_class":"geometry",
        "evidence_class":"uncalibrated_model_output"})
    attrs = {
        "model_version":"v0.6.0-development", "axis_order":"time,z,y,x",
        "flatten_order":"C", "geometry_hash":grid.geometry_hash,
        "config_hash":config_hash(result.config), "grid_shape":json.dumps(grid.shape),
        "random_seed":result.config.simulation.random_seed,
        "linear_solver":result.config.solver_3d.linear_solver,
        "relative_tolerance":result.config.solver_3d.relative_tolerance,
        "absolute_tolerance":result.config.solver_3d.absolute_tolerance,
        "evidence_label":"Uncalibrated 3D model output; Not experimental data",
        "team_wet_lab_rows":0, "public_calibration_status":"not calibrated",
    }
    return xr.Dataset(variables, coords=coordinates, attrs=attrs)


def _git(root: Path, *arguments: str) -> str:
    try:
        return subprocess.check_output(["git", *arguments], cwd=str(root), text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def require_gate_d(validation_report: Path) -> Dict[str, Any]:
    path = Path(validation_report)
    if not path.exists():
        raise RuntimeError("Formal 3D output is locked until Gate D validation exists")
    report = json.loads(path.read_text(encoding="utf-8"))
    if not report.get("gate_d_passed"):
        raise RuntimeError("Formal 3D output is locked because Gate D did not pass")
    return report


def save_formal_run_3d(result, base: Path, validation_report: Path,
                       run_id: str = None) -> Path:
    """Atomically write the required Zarr run directory; never falls back to CSV."""
    require_gate_d(validation_report)
    xr, Blosc = _optional_storage_modules()
    run_id = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    final = Path(base) / run_id
    if final.exists():
        raise FileExistsError("Refusing to overwrite existing 3D run: {}".format(final))
    temporary = final.parent / ("." + run_id + ".tmp-" + uuid.uuid4().hex)
    temporary.mkdir(parents=True, exist_ok=False)
    try:
        dataset = result_to_dataset(result)
        compressor = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
        encoding = {name:{"compressor":compressor} for name in dataset.data_vars}
        dataset.to_zarr(str(temporary / "fields.zarr"), mode="w", encoding=encoding,
                        consolidated=True)
        reopened = xr.open_zarr(str(temporary / "fields.zarr"), consolidated=True)
        if dict(reopened.sizes) != dict(dataset.sizes):
            raise RuntimeError("Zarr verification failed: dimension mismatch")
        reopened.close()
        root = Path(__file__).resolve().parents[1]
        commit = _git(root, "rev-parse", "HEAD")
        dirty = bool(_git(root, "status", "--short"))
        payloads = {
            "summary.json":result.summary, "diagnostics.json":result.diagnostics,
            "performance.json":result.performance, "config.json":result.config.to_dict(),
            "geometry.json":{"shape":result.geometry.shape,"axis_order":["z","y","x"],
                             "geometry_hash":result.geometry.geometry_hash},
            "boundary_conditions.json":result.config.boundary_3d.__dict__,
            "run_manifest.json":{"model_version":"v0.6.0-development","git_commit":commit,
                "dirty_worktree":dirty,"config_hash":config_hash(result.config),
                "geometry_hash":result.geometry.geometry_hash,"grid_shape":result.geometry.shape,
                "axis_order":["z","y","x"],"random_seed":result.config.simulation.random_seed,
                "solver":result.config.solver_3d.linear_solver,
                "evidence_label":"Uncalibrated 3D model output; Not experimental data",
                "public_calibration_status":"not calibrated","gate_d_passed":True},
        }
        for name, payload in payloads.items():
            (temporary/name).write_text(json.dumps(payload,indent=2),encoding="utf-8")
        (temporary/"checkpoints").mkdir(); (temporary/"figures").mkdir()
        final.parent.mkdir(parents=True, exist_ok=True)
        temporary.replace(final)
    except Exception:
        shutil.rmtree(str(temporary), ignore_errors=True)
        raise
    return final

"""Surface-resolved solid partition and fixed-grid aperture evolution."""

from dataclasses import dataclass
import numpy as np
from .config import ModelConfig
from .grid_3d import StructuredGrid3D
from .state import S

@dataclass(frozen=True)
class DepositionResult3D:
    aperture_m: np.ndarray
    closure: np.ndarray
    sealed_mask: np.ndarray
    calcite_wall_volume_m3: np.ndarray
    csh_wall_volume_m3: np.ndarray
    wall_solid_volume_m3: np.ndarray
    bulk_solid_volume_m3: np.ndarray
    total_solid_volume_m3: np.ndarray
    area_weighted_closure: float
    open_volume_closure: float
    sealed_area_fraction: float

def update_aperture_3d(state, grid: StructuredGrid3D, config: ModelConfig,
                       initial_aperture_m=None, minimum_aperture_m=0.0,
                       closure_threshold_m=1e-9) -> DepositionResult3D:
    values=np.asarray(state).reshape(grid.shape+(state.shape[-1],))
    b0=np.full((grid.shape[0],grid.shape[2]),config.transport.crack_width_mm*1e-3) if initial_aperture_m is None else np.asarray(initial_aperture_m,float)
    calcite_fraction=values[...,S['calcite_mol_m3']]*config.chemistry.calcite_molar_mass_kg_mol/config.chemistry.calcite_density_kg_m3
    csh_fraction=values[...,S['csh_volume_fraction']]
    calcite_volume=np.sum(calcite_fraction*grid.cell_volume_m3,axis=1)
    csh_volume=np.sum(csh_fraction*grid.cell_volume_m3,axis=1)
    total=calcite_volume+csh_volume; fraction=config.chemistry.wall_deposition_fraction
    calcite_wall=calcite_volume*fraction; csh_wall=csh_volume*fraction
    wall=calcite_wall+csh_wall; bulk=total-wall
    projected_area=np.sum(grid.face_area_y_m2,axis=1)/grid.shape[1]
    aperture=np.maximum(b0-wall/np.maximum(projected_area,1e-300),minimum_aperture_m)
    closure=np.clip(1-aperture/np.maximum(b0,1e-300),0,1)
    weights=projected_area
    area_closure=1-float(np.sum(aperture*weights)/np.sum(b0*weights))
    volume_closure=1-float(np.sum(aperture*weights)/np.sum(b0*weights))
    sealed=aperture<=closure_threshold_m
    return DepositionResult3D(aperture,closure,sealed,calcite_wall,csh_wall,wall,bulk,total,
                              area_closure,volume_closure,float(np.sum(weights*sealed)/np.sum(weights)))

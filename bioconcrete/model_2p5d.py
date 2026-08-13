"""Aperture-averaged production geometry and topology-safe flow metrics."""
from dataclasses import dataclass
import numpy as np
from .config import ModelConfig

@dataclass(frozen=True)
class FlowResult2p5D:
    pressure_pa: np.ndarray
    flux_x_m2_s: np.ndarray
    total_flow_m3_s: float
    relative_transmissivity: float
    mean_aperture_cubic_proxy: float
    is_through_flow: bool

def solve_flow_2p5d(aperture_m: np.ndarray, config: ModelConfig,
                    inlet_pressure_pa: float=1.0, outlet_pressure_pa: float=0.0,
                    viscosity_pa_s: float=1e-3, reference_aperture_m=None) -> FlowResult2p5D:
    b=np.asarray(aperture_m,float)
    if b.ndim!=2 or np.any(b<0): raise ValueError('aperture must be a nonnegative (z,x) field')
    b0=config.transport.crack_width_mm*1e-3 if reference_aperture_m is None else float(reference_aperture_m)
    proxy=float(np.mean((b/max(b0,1e-300))**3))
    if config.geometry_3d.topology=='blind_crack':
        return FlowResult2p5D(np.full_like(b,np.nan),np.zeros_like(b),float('nan'),float('nan'),proxy,False)
    if viscosity_pa_s<=0: raise ValueError('viscosity must be positive')
    nx=b.shape[1]; length=config.transport.crack_length_mm*1e-3; dx=length/nx
    pressure=np.broadcast_to(np.linspace(inlet_pressure_pa,outlet_pressure_pa,nx),b.shape).copy()
    gradient=(outlet_pressure_pa-inlet_pressure_pa)/length
    flux=-b**3/(12*viscosity_pa_s)*gradient
    dz=config.transport.crack_depth_mm*1e-3/b.shape[0]
    total=float(np.sum(flux[:,-1]*dz))
    reference=float((b0**3/(12*viscosity_pa_s))*(-gradient)*config.transport.crack_depth_mm*1e-3)
    relative=total/reference if reference else float('nan')
    return FlowResult2p5D(pressure,flux,total,relative,proxy,True)

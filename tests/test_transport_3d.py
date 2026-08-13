import unittest
import numpy as np
from bioconcrete.boundary_3d import BoundaryCondition3D
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import rectangular_grid_3d
from bioconcrete.transport_3d import transport_step_3d

class Transport3DTests(unittest.TestCase):
    def grid(self,n=7):
        c=ModelConfig(); c.geometry_3d.nx=n; c.geometry_3d.ny=3; c.geometry_3d.nz=5
        return rectangular_grid_3d(c)
    def test_constant_no_flux_and_zero_gradient(self):
        g=self.grid(); field=np.full(g.shape,2.0)
        result,d=transport_step_3d(field,g,10,1e-9)
        np.testing.assert_allclose(result,field,atol=1e-12); self.assertAlmostEqual(d.boundary_rate_before,0,places=18)
    def test_real_z_gradient_diffuses(self):
        g=self.grid(); field=np.zeros(g.shape); field[0]=1
        result,_=transport_step_3d(field,g,10000,1e-9)
        self.assertGreater(result[1].mean(),0); self.assertLess(result[0].mean(),1)
        self.assertGreater(np.ptp(result.mean(axis=(1,2))),0)
    def test_z_uniform_stays_uniform(self):
        g=self.grid(); field=np.broadcast_to(np.linspace(0,1,g.shape[2]),g.shape).copy()
        result,_=transport_step_3d(field,g,10,1e-9)
        np.testing.assert_allclose(result[0],result[-1],atol=1e-12)
    def test_boundary_inventory_change_matches_implicit_flux(self):
        g=self.grid(); field=np.zeros(g.shape); bc={'x_min':BoundaryCondition3D('dirichlet',1)}
        dt=50; result,d=transport_step_3d(field,g,dt,1e-9,bc)
        change=np.sum((result-field)*g.cell_volume_m3)
        self.assertAlmostEqual(change,dt*d.boundary_rate_after,delta=max(abs(change)*1e-9,1e-20))
        self.assertAlmostEqual(sum(d.boundary_rates_after.values()),d.boundary_rate_after,
                               delta=max(abs(d.boundary_rate_after)*1e-9,1e-20))
    def test_direct_and_iterative_agree(self):
        g=self.grid(); rng=np.random.RandomState(2); field=rng.rand(*g.shape)
        direct,_=transport_step_3d(field,g,10,1e-9,linear_solver='direct')
        iterative,_=transport_step_3d(field,g,10,1e-9,linear_solver='cg',relative_tolerance=1e-11,absolute_tolerance=1e-13)
        np.testing.assert_allclose(direct,iterative,rtol=1e-8,atol=1e-10)

if __name__ == '__main__': unittest.main()

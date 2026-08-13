import math
import unittest
import numpy as np
from bioconcrete.boundary_3d import BoundaryCondition3D
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import rectangular_grid_3d
from bioconcrete.transport_3d import build_transport_operator_3d

class ManufacturedSolution3DTests(unittest.TestCase):
    def error(self,n):
        c=ModelConfig(); c.transport.crack_length_mm=10; c.transport.crack_width_mm=10; c.transport.crack_depth_mm=10
        c.geometry_3d.nx=n; c.geometry_3d.ny=n; c.geometry_3d.nz=n
        g=rectangular_grid_3d(c); zz,yy,xx=g.cell_coordinates(); length=.01
        exact=np.sin(math.pi*xx/length)*np.sin(math.pi*yy/length)*np.sin(math.pi*zz/length)
        laplace=-3*(math.pi/length)**2*exact
        bc={face:BoundaryCondition3D('dirichlet',0) for face in ('x_min','x_max','y_min','y_max','z_min','z_max')}
        operator,source=build_transport_operator_3d(g,1.0,bc)
        residual=g.reshape(operator@g.flatten(exact)+source)-laplace
        return math.sqrt(np.mean(residual**2))/math.sqrt(np.mean(laplace**2))
    def test_diffusion_spatial_order_is_second(self):
        coarse=self.error(8); fine=self.error(16)
        order=math.log(coarse/fine,2)
        self.assertGreaterEqual(order,1.8)

if __name__ == '__main__': unittest.main()

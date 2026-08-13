import unittest
import numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import rectangular_grid_3d

class Grid3DTests(unittest.TestCase):
    def test_volume_area_order_and_round_trip(self):
        c=ModelConfig(); c.geometry_3d.nx=7; c.geometry_3d.ny=3; c.geometry_3d.nz=5
        g=rectangular_grid_3d(c)
        self.assertEqual(g.shape,(5,3,7)); self.assertEqual(g.size,105)
        self.assertAlmostEqual(np.sum(g.cell_volume_m3), .1*.0003*.02)
        self.assertAlmostEqual(np.sum(g.face_area_y_m2[:,0,:]), .1*.02)
        a=np.arange(g.size).reshape(g.shape); np.testing.assert_array_equal(g.reshape(g.flatten(a)),a)
        self.assertEqual(len(g.geometry_hash),64)

if __name__ == '__main__': unittest.main()

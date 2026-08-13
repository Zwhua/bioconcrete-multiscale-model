import unittest
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import rectangular_grid_3d

class Geometry3DTests(unittest.TestCase):
    def test_legacy_config_uses_defaults(self):
        c=ModelConfig.load(); self.assertEqual(c.geometry_3d.topology,'blind_crack')
    def test_cell_centres_stay_inside_physical_domain(self):
        c=ModelConfig(); g=rectangular_grid_3d(c)
        self.assertGreater(g.x_m.min(),0); self.assertLess(g.x_m.max(),.1)
        self.assertLess(g.y_m.max(),.0003); self.assertLess(g.z_m.max(),.02)

if __name__ == '__main__': unittest.main()

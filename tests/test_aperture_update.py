import unittest,numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.deposition_3d import update_aperture_3d
from bioconcrete.grid_3d import rectangular_grid_3d
from bioconcrete.state import STATE_NAMES
class ApertureTests(unittest.TestCase):
 def test_zero_deposition_keeps_aperture(self):
  c=ModelConfig();g=rectangular_grid_3d(c);s=np.zeros((g.size,len(STATE_NAMES)));r=update_aperture_3d(s,g,c)
  np.testing.assert_allclose(r.aperture_m,c.transport.crack_width_mm*1e-3)
 def test_zero_wall_fraction_prevents_closure(self):
  c=ModelConfig();c.chemistry.wall_deposition_fraction=0;g=rectangular_grid_3d(c);s=np.ones((g.size,len(STATE_NAMES)));r=update_aperture_3d(s,g,c)
  self.assertEqual(r.area_weighted_closure,0)
if __name__=='__main__':unittest.main()

import unittest, numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.deposition_3d import update_aperture_3d
from bioconcrete.grid_3d import rectangular_grid_3d
from bioconcrete.state import S,STATE_NAMES
class Deposition3DTests(unittest.TestCase):
 def test_partition_and_separate_contributions(self):
  c=ModelConfig(); c.geometry_3d.nx=3;c.geometry_3d.ny=3;c.geometry_3d.nz=2;g=rectangular_grid_3d(c)
  s=np.zeros((g.size,len(STATE_NAMES)));s[:,S['calcite_mol_m3']]=2;s[:,S['csh_volume_fraction']]=.001
  r=update_aperture_3d(s,g,c)
  np.testing.assert_allclose(r.wall_solid_volume_m3+r.bulk_solid_volume_m3,r.total_solid_volume_m3)
  np.testing.assert_allclose(r.calcite_wall_volume_m3+r.csh_wall_volume_m3,r.wall_solid_volume_m3)
  self.assertTrue(np.all((r.closure>=0)&(r.closure<=1)))
if __name__=='__main__':unittest.main()

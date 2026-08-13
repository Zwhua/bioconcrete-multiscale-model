import unittest,numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.model_3d import simulate_3d
class Model3DTests(unittest.TestCase):
 def test_minimal_model_runs_and_has_true_z_difference(self):
  c=ModelConfig();c.geometry_3d.nx=9;c.geometry_3d.ny=3;c.geometry_3d.nz=5;c.simulation.days=.0007;c.simulation.reaction_step_h=.02;c.output_3d.save_every_days=.0007
  c.geometry_3d.capsule_depth_mode='surface';r=simulate_3d(c)
  self.assertEqual(r.state.shape[1:4],(5,3,9));self.assertTrue(r.diagnostics['finite']);self.assertTrue(r.diagnostics['nonnegative'])
  self.assertGreater(r.diagnostics['z_oxygen_range_mol_m3'],0)
  self.assertIn('not experimental data',r.diagnostics['evidence_label'])
  self.assertEqual(r.aperture_history_m.shape,(2,5,9))
if __name__=='__main__':unittest.main()

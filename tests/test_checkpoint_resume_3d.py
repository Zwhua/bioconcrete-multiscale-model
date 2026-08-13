import tempfile,unittest
from pathlib import Path
import numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.model_3d import simulate_3d

class CheckpointResume3DTests(unittest.TestCase):
 def test_interrupted_and_continuous_results_match(self):
  c=ModelConfig();c.geometry_3d.nx=3;c.geometry_3d.ny=2;c.geometry_3d.nz=2
  c.simulation.days=.0004;c.simulation.reaction_step_h=.0048;c.output_3d.save_every_days=.0002
  c.solver_3d.checkpoint_interval_steps=1
  continuous=simulate_3d(c)
  with tempfile.TemporaryDirectory() as d:
   checkpoint=Path(d)/'state.npz'
   simulate_3d(c,checkpoint=checkpoint,stop_after_s=c.simulation.days*86400/2)
   resumed=simulate_3d(c,checkpoint=checkpoint,resume=True)
  np.testing.assert_allclose(resumed.state[-1],continuous.state[-1],rtol=1e-10,atol=1e-12)
  self.assertEqual(resumed.diagnostics['conservation'],continuous.diagnostics['conservation'])

if __name__=='__main__':unittest.main()

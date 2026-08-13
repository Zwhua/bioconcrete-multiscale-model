import tempfile,unittest,numpy as np
from pathlib import Path
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import rectangular_grid_3d
from bioconcrete.io_3d import read_checkpoint,write_checkpoint
from bioconcrete.state import STATE_NAMES
class Checkpoint3DTests(unittest.TestCase):
 def test_roundtrip_and_hash_rejection(self):
  c=ModelConfig();c.geometry_3d.nx=3;c.geometry_3d.ny=2;c.geometry_3d.nz=2;g=rectangular_grid_3d(c);s=np.arange(g.size*len(STATE_NAMES)).reshape(g.size,-1)
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'state.npz';write_checkpoint(p,12,s,c,g);t,r=read_checkpoint(p,c,g);self.assertEqual(t,12);np.testing.assert_array_equal(r,s);self.assertFalse(Path(str(p)+'.tmp').exists())
   write_checkpoint(p,13,s,c,g,{'carbon':1},4,2,1);t,r,m=read_checkpoint(p,c,g,full=True)
   self.assertEqual(m['ledger'],{'carbon':1});self.assertEqual(m['step'],4);self.assertEqual(m['output_index'],2)
   c.simulation.days+=1
   with self.assertRaises(ValueError):read_checkpoint(p,c,g)
if __name__=='__main__':unittest.main()

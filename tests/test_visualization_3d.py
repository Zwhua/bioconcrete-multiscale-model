import tempfile,unittest
from pathlib import Path
from bioconcrete.visualization_3d import render_formal_3d
class Visualization3DTests(unittest.TestCase):
 def test_gate_lock_prevents_formal_figures(self):
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaises(RuntimeError):render_formal_3d(Path(d),Path(d)/'missing.json')
if __name__=='__main__':unittest.main()

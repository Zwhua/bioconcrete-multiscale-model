import tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from bioconcrete.io_3d import require_gate_d,result_to_dataset

class IO3DTests(unittest.TestCase):
 def test_gate_lock(self):
  with tempfile.TemporaryDirectory() as d:
   with self.assertRaises(RuntimeError):require_gate_d(Path(d)/'missing.json')
 def test_missing_optional_dependency_is_explicit(self):
  with patch.dict('sys.modules',{'xarray':None}):
   with self.assertRaisesRegex(RuntimeError,"three-d"):
    from bioconcrete.io_3d import _optional_storage_modules;_optional_storage_modules()
if __name__=='__main__':unittest.main()

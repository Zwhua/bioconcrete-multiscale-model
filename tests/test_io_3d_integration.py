import tempfile,unittest
from pathlib import Path
from bioconcrete.config import ModelConfig
from bioconcrete.model_3d import simulate_3d
from bioconcrete.io_3d import result_to_dataset

class IO3DIntegrationTests(unittest.TestCase):
 def test_dataset_dimensions_and_metadata(self):
  try: import xarray,zarr,numcodecs
  except ImportError: self.skipTest('three-d extra is not installed')
  c=ModelConfig();c.geometry_3d.nx=3;c.geometry_3d.ny=2;c.geometry_3d.nz=2;c.simulation.days=.0001;c.output_3d.save_every_days=.0001;c.simulation.reaction_step_h=.0024
  result=simulate_3d(c);ds=result_to_dataset(result)
  self.assertEqual(ds['oxygen_mol_m3'].dims,('time','z','y','x'));self.assertEqual(ds['aperture_m'].dims,('time','z','x'))
  self.assertIn('evidence_class',ds['ph'].attrs);self.assertIn('Not experimental data',ds.attrs['evidence_label'])
if __name__=='__main__':unittest.main()

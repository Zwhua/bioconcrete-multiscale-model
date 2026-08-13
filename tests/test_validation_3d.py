import json,tempfile,unittest
from pathlib import Path
from bioconcrete.validation_3d import run_validation_3d,validation_config

class Validation3DTests(unittest.TestCase):
    def test_smoke_never_claims_gate_d(self):
        with tempfile.TemporaryDirectory() as d:
            result=run_validation_3d(Path(d),full=False)
            self.assertFalse(result['gate_d_passed'])
            self.assertTrue(result['checks']['conservation_passed'])
            self.assertTrue(result['checks']['reduction_passed'])
            self.assertTrue((Path(d)/'convergence.csv').exists())
    def test_frozen_scenario_is_nonzero_and_documented(self):
        c=validation_config();self.assertEqual(c.simulation.days,1);self.assertTrue(c.simulation.closed_system)

if __name__=='__main__': unittest.main()

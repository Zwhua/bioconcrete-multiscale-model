import unittest,numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.model_2p5d import solve_flow_2p5d
class Flow2p5DTests(unittest.TestCase):
 def test_blind_never_reports_fake_flow(self):
  c=ModelConfig();r=solve_flow_2p5d(np.full((3,5),3e-4),c);self.assertFalse(r.is_through_flow);self.assertTrue(np.isnan(r.total_flow_m3_s))
 def test_uniform_through_flow_and_zero_pressure(self):
  c=ModelConfig();c.geometry_3d.topology='through_crack';b=np.full((3,5),3e-4)
  r=solve_flow_2p5d(b,c);self.assertAlmostEqual(r.relative_transmissivity,1);self.assertGreater(r.total_flow_m3_s,0)
  z=solve_flow_2p5d(b,c,0,0);self.assertEqual(z.total_flow_m3_s,0)
 def test_spatial_distribution_changes_flow(self):
  c=ModelConfig();c.geometry_3d.topology='through_crack';a=np.full((2,4),2e-4);b=a.copy();b[0]=1e-4;b[1]=3e-4
  qa=solve_flow_2p5d(a,c).total_flow_m3_s;qb=solve_flow_2p5d(b,c).total_flow_m3_s
  self.assertGreater(abs(qb-qa)/qa,0.1)
if __name__=='__main__':unittest.main()

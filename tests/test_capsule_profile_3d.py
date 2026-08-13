import unittest
import numpy as np
from bioconcrete.config import ModelConfig
from bioconcrete.grid_3d import capsule_profile_3d, rectangular_grid_3d

class CapsuleProfile3DTests(unittest.TestCase):
    def inventory(self,c):
        g=rectangular_grid_3d(c); p,centres=capsule_profile_3d(c,g)
        return np.sum(p*g.cell_volume_m3), centres
    def test_seed_and_volume_normalisation(self):
        c=ModelConfig(); a,ca=self.inventory(c); b,cb=self.inventory(c)
        self.assertAlmostEqual(a,.1*.0003*.02); np.testing.assert_array_equal(ca,cb)
    def test_refinement_and_width_preserve_fixed_inventory_factor(self):
        c=ModelConfig(); a,_=self.inventory(c)
        c.geometry_3d.nx*=2; c.geometry_3d.nz*=2; b,_=self.inventory(c)
        self.assertAlmostEqual(a,b)
        c.transport.crack_width_mm=.15; narrow,_=self.inventory(c)
        self.assertAlmostEqual(narrow,a/2)

if __name__ == '__main__': unittest.main()

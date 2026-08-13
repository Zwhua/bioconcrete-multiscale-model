import unittest
from bioconcrete.config import ModelConfig
from bioconcrete.model import _initial_state
from bioconcrete.state import S
import numpy as np

class InitialCarbon3DTests(unittest.TestCase):
    def test_explicit_initial_dic_is_loaded_and_hydrated(self):
        c=ModelConfig();c.environment.inorganic_carbon_initial_mol_m3=.15
        state=_initial_state(c,np.ones(2),'3d')
        np.testing.assert_allclose(state[:,S['inorganic_carbon_mol_m3']],.15)
        self.assertTrue(np.all(state[:,S['hydrated_carbon_mol_m3']]>0))

if __name__=='__main__':unittest.main()

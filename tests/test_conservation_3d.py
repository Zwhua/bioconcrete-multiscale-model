import unittest
from bioconcrete.config import ModelConfig
from bioconcrete.model_3d import simulate_3d

def small_closed():
    c=ModelConfig();c.geometry_3d.nx=5;c.geometry_3d.ny=3;c.geometry_3d.nz=3
    c.simulation.days=.002;c.simulation.reaction_step_h=.024;c.output_3d.save_every_days=.002
    c.simulation.closed_system=True;c.environment.exposure='continuous'
    c.environment.oxygen_initial_mol_m3=.25
    return c

class Conservation3DTests(unittest.TestCase):
    def test_closed_carbon_and_calcium(self):
        r=simulate_3d(small_closed())
        self.assertLess(r.diagnostics['conservation']['carbon_mol']['relative_error'],.005)
        self.assertLess(r.diagnostics['conservation']['calcium_mol']['relative_error'],.005)
    def test_open_balance_records_boundary_integrals(self):
        c=small_closed();c.simulation.closed_system=False;c.environment.oxygen_initial_mol_m3=0
        r=simulate_3d(c)
        self.assertIn('boundary_integral',r.diagnostics['conservation']['carbon_mol'])
        self.assertLess(r.diagnostics['conservation']['carbon_mol']['relative_error'],.005)

if __name__=='__main__': unittest.main()

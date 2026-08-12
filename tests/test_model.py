import copy
import unittest

import numpy as np

from bioconcrete.chemistry import carbonate_fractions
from bioconcrete.config import ModelConfig
from bioconcrete.model import simulate_0d, simulate_1d, simulate_2d


class ModelTests(unittest.TestCase):
    def short_config(self):
        config = ModelConfig()
        config.simulation.days = 0.05
        config.simulation.output_interval_days = 0.025
        config.simulation.reaction_step_h = 1.2
        return config

    def test_carbonate_fractions_sum_to_one(self):
        ph = np.linspace(6.0, 13.0, 30)
        fractions = carbonate_fractions(ph, 30.0)
        np.testing.assert_allclose(sum(fractions), np.ones_like(ph), rtol=1e-12, atol=1e-12)

    def test_closed_zero_dimensional_balance_and_ammonia_invariant(self):
        config = self.short_config()
        config.simulation.closed_system = True
        result = simulate_0d(config)
        changes = result.diagnostics["relative_balance_change"]
        self.assertLess(changes["carbon"], 5e-3)
        self.assertLess(changes["calcium"], 5e-3)
        self.assertTrue(result.diagnostics["nonnegative"])
        self.assertTrue(result.diagnostics["ammonia_free"])
        self.assertEqual(result.summary["ammonium_mol_m3_max"], 0.0)

    def test_no_active_material_means_no_mineralization(self):
        config = self.short_config()
        config.simulation.closed_system = True
        config.kinetics.capsule_calcium_lactate_mol_m3 = 0.0
        config.kinetics.spore_density_rel = 0.0
        config.kinetics.active_density_rel = 0.0
        result = simulate_0d(config)
        self.assertLess(result.summary["calcite_mol_m3_mean"], 1e-12)
        self.assertLess(result.summary["mean_healing_ratio"], 1e-12)

    def test_one_dimensional_solver_has_spatial_coordinates(self):
        config = self.short_config()
        config.transport.nx_1d = 9
        config.transport.capsule_count_1d = 2
        result = simulate_1d(config)
        self.assertIn("x_mm", result.frame.columns)
        self.assertEqual(result.diagnostics["grid_shape"], [9])
        self.assertTrue(result.diagnostics["true_spatial_solver"])
        self.assertTrue((result.frame.select_dtypes(include=[np.number]) >= -1e-10).all().all())

    def test_two_dimensional_solver_is_not_interpolated_one_dimensional_output(self):
        config = self.short_config()
        config.transport.nx_2d = 7
        config.transport.ny_2d = 5
        config.transport.capsule_count_2d = 4
        result = simulate_2d(config)
        self.assertIn("x_mm", result.frame.columns)
        self.assertIn("y_mm", result.frame.columns)
        self.assertEqual(result.diagnostics["grid_shape"], [5, 7])
        self.assertEqual(len(result.frame[result.frame["time_d"] == result.frame["time_d"].max()]), 35)
        self.assertTrue(result.diagnostics["true_spatial_solver"])


if __name__ == "__main__":
    unittest.main()

import copy
import unittest

import numpy as np

from bioconcrete.config import ModelConfig
from bioconcrete.model import (
    S, STATE_NAMES, _initial_state, _inventory, repair_metrics,
    simulate_0d, simulate_1d, simulate_2d,
)


class PhysicalLimitTests(unittest.TestCase):
    def short_config(self):
        config = ModelConfig()
        config.simulation.days = 0.05
        config.simulation.output_interval_days = 0.025
        config.simulation.reaction_step_h = 1.2
        return config

    def test_zero_dosage_zeros_all_repair_agent_inventories(self):
        config = self.short_config()
        config.kinetics.agent_dosage_multiplier = 0.0
        state = _initial_state(config, np.ones(1), "0d")
        inventory = _inventory(state, config, "0d")
        self.assertEqual(inventory["capsule_calcium_lactate_mol"], 0.0)
        self.assertEqual(inventory["spore_inventory_rel_m3"], 0.0)
        self.assertEqual(inventory["active_inventory_rel_m3"], 0.0)
        self.assertEqual(inventory["remaining_csh_payload_m3"], 0.0)

    def test_double_dosage_doubles_complete_agent_inventory(self):
        base = self.short_config()
        doubled = copy.deepcopy(base)
        doubled.kinetics.agent_dosage_multiplier = 2.0
        first = _inventory(_initial_state(base, np.ones(1), "0d"), base, "0d")
        second = _inventory(_initial_state(doubled, np.ones(1), "0d"), doubled, "0d")
        for name in (
            "capsule_calcium_lactate_mol", "spore_inventory_rel_m3",
            "active_inventory_rel_m3", "remaining_csh_payload_m3",
        ):
            self.assertAlmostEqual(second[name], 2.0 * first[name], places=15)

    def test_fixed_total_inventory_is_independent_of_crack_width(self):
        narrow = self.short_config()
        wide = copy.deepcopy(narrow)
        wide.transport.crack_width_mm = 0.5
        narrow.transport.crack_width_mm = 0.1
        narrow_inventory = _inventory(
            _initial_state(narrow, np.ones(1), "0d"), narrow, "0d"
        )["capsule_calcium_lactate_mol"]
        wide_inventory = _inventory(
            _initial_state(wide, np.ones(1), "0d"), wide, "0d"
        )["capsule_calcium_lactate_mol"]
        self.assertAlmostEqual(narrow_inventory, wide_inventory, places=15)

    def test_zero_spores_with_csh_produces_only_nonbiological_fill(self):
        config = self.short_config()
        config.kinetics.spore_density_rel = 0.0
        config.kinetics.active_density_rel = 0.0
        config.chemistry.calcite_rate_mol_m3_s = 0.0
        result = simulate_0d(config)
        self.assertGreater(result.frame["csh_solid_volume_m3"].iloc[-1], 0.0)
        self.assertEqual(result.summary["mean_calcite_closure_contribution"], 0.0)
        self.assertGreater(result.summary["mean_csh_closure_contribution"], 0.0)

    def test_zero_calcium_sources_prevent_calcite(self):
        config = self.short_config()
        config.kinetics.capsule_calcium_lactate_mol_m3 = 0.0
        config.chemistry.portlandite_mol_m3 = 0.0
        config.kinetics.capsule_csh_volume_fraction = 0.0
        result = simulate_0d(config)
        self.assertLess(result.summary["calcite_mol_m3_mean"], 1e-12)

    def test_zero_carbon_inventory_prevents_calcite_in_closed_system(self):
        config = self.short_config()
        config.simulation.closed_system = True
        config.kinetics.agent_dosage_multiplier = 0.0
        config.kinetics.capsule_csh_volume_fraction = 0.0
        result = simulate_0d(config)
        self.assertLess(result.summary["calcite_mol_m3_mean"], 1e-12)

    def test_wide_crack_closes_less_for_same_inventory(self):
        narrow = self.short_config()
        wide = copy.deepcopy(narrow)
        narrow.transport.crack_width_mm = 0.1
        wide.transport.crack_width_mm = 0.5
        self.assertGreater(
            simulate_0d(narrow).summary["mean_crack_closure_ratio"],
            simulate_0d(wide).summary["mean_crack_closure_ratio"],
        )

    def test_zero_wall_fraction_has_no_wall_volume_or_closure(self):
        config = self.short_config()
        config.chemistry.wall_deposition_fraction = 0.0
        result = simulate_0d(config)
        self.assertEqual(result.frame["wall_solid_volume_m3"].max(), 0.0)
        self.assertEqual(result.summary["mean_crack_closure_ratio"], 0.0)

    def test_outputs_are_finite_nonnegative_and_bounded(self):
        result = simulate_0d(self.short_config())
        numeric = result.frame.select_dtypes(include=[np.number]).to_numpy()
        self.assertTrue(np.isfinite(numeric).all())
        concentration_columns = [name for name in STATE_NAMES if name != "tracked_ph"]
        self.assertTrue((result.frame[concentration_columns] >= -1e-12).all().all())
        self.assertTrue(result.frame["crack_closure_ratio"].between(0.0, 1.0).all())

    def test_seeded_two_dimensional_result_is_reproducible(self):
        config = self.short_config()
        config.transport.nx_2d = 7
        config.transport.ny_2d = 3
        first = simulate_2d(config).frame
        second = simulate_2d(copy.deepcopy(config)).frame
        np.testing.assert_allclose(
            first["crack_closure_ratio"], second["crack_closure_ratio"],
            rtol=0.0, atol=0.0,
        )

    def test_grid_refinement_preserves_initial_total_inventory(self):
        coarse = self.short_config()
        fine = copy.deepcopy(coarse)
        coarse.transport.nx_1d = 9
        fine.transport.nx_1d = 17
        coarse_result = simulate_1d(coarse)
        fine_result = simulate_1d(fine)
        self.assertAlmostEqual(
            coarse_result.diagnostics["initial_inventory"]["capsule_calcium_lactate_mol"],
            fine_result.diagnostics["initial_inventory"]["capsule_calcium_lactate_mol"],
            places=15,
        )

    def test_calcite_and_csh_contributions_sum_before_saturation(self):
        config = self.short_config()
        state = np.zeros((1, len(STATE_NAMES)))
        state[0, S["calcite_mol_m3"]] = 1.0
        state[0, S["csh_volume_fraction"]] = 1e-4
        metrics = repair_metrics(state, config)
        self.assertAlmostEqual(
            metrics["calcite_closure_contribution"][0]
            + metrics["csh_closure_contribution"][0],
            metrics["crack_closure_ratio"][0],
            places=12,
        )


if __name__ == "__main__":
    unittest.main()

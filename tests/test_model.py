import copy
import unittest

import numpy as np

from bioconcrete.chemistry import alkalinity_from_ph, carbonate_fractions, ph_from_alkalinity
from bioconcrete.config import ModelConfig
from bioconcrete.model import (
    S, STATE_NAMES, _capsule_profile_1d, _capsule_profile_2d, _initial_state,
    _inventory, repair_metrics, simulate_0d, simulate_1d, simulate_2d,
)


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

    def test_dynamic_ph_charge_balance_round_trip(self):
        carbon = np.array([0.0, 10.0, 100.0])
        expected = np.array([9.0, 10.5, 12.0])
        alkalinity = alkalinity_from_ph(carbon, expected, 30.0)
        actual = ph_from_alkalinity(carbon, alkalinity, 30.0, 6.0, 13.5)
        np.testing.assert_allclose(actual, expected, atol=1e-10)

    def test_dynamic_and_fixed_ph_modes_are_explicit(self):
        dynamic = self.short_config()
        fixed = copy.deepcopy(dynamic)
        fixed.simulation.ph_mode = "fixed"
        dynamic_result = simulate_0d(dynamic)
        fixed_result = simulate_0d(fixed)
        self.assertTrue((dynamic_result.frame["ph"] >= dynamic.environment.ph_minimum).all())
        np.testing.assert_allclose(fixed_result.frame["ph"], fixed.environment.ph)

    def test_geometry_outputs_and_thickness_mass_scaling(self):
        thin = self.short_config()
        thick = copy.deepcopy(thin)
        thick.transport.crack_depth_mm *= 2.0
        thin_result = simulate_0d(thin)
        thick_result = simulate_0d(thick)
        self.assertIn("solid_fill_fraction", thin_result.frame)
        self.assertIn("wall_deposition_thickness_mm", thin_result.frame)
        self.assertIn("crack_closure_ratio", thin_result.frame)
        self.assertAlmostEqual(
            thick_result.frame["cell_volume_m3"].iloc[-1] / thin_result.frame["cell_volume_m3"].iloc[-1],
            2.0, places=7,
        )
        self.assertAlmostEqual(
            thick_result.summary["calcite_mol_m3_mean"], thin_result.summary["calcite_mol_m3_mean"], places=7
        )
        self.assertAlmostEqual(
            thick_result.summary["mean_crack_closure_ratio"],
            thin_result.summary["mean_crack_closure_ratio"], places=7,
        )

    def test_wall_area_controls_deposition_thickness(self):
        config = self.short_config()
        state = np.zeros((1, len(STATE_NAMES)))
        state[0, S["calcite_mol_m3"]] = 100.0
        small_area = repair_metrics(state, config, np.array([1.0e-9]), np.array([2.0e-6]))
        large_area = repair_metrics(state, config, np.array([1.0e-9]), np.array([4.0e-6]))
        self.assertAlmostEqual(
            small_area["wall_deposition_thickness_mm"][0] /
            large_area["wall_deposition_thickness_mm"][0],
            2.0,
        )

    def test_wall_and_nonwall_solid_volumes_partition_total(self):
        config = self.short_config()
        state = np.zeros((1, len(STATE_NAMES)))
        state[0, S["calcite_mol_m3"]] = 100.0
        metrics = repair_metrics(state, config, np.array([1.0e-9]), np.array([2.0e-6]))
        np.testing.assert_allclose(
            metrics["wall_solid_volume_m3"] + metrics["nonwall_solid_volume_m3"],
            metrics["total_solid_volume_m3"],
        )

    def test_initial_capsule_inventory_matches_across_scales(self):
        config = self.short_config()
        config.transport.out_of_plane_thickness_mm = config.transport.crack_depth_mm
        _, profile_1d = _capsule_profile_1d(config)
        _, _, profile_2d = _capsule_profile_2d(config)
        states = {
            "0d": _initial_state(config, np.ones(1)),
            "1d": _initial_state(config, profile_1d),
            "2d": _initial_state(config, profile_2d.ravel()),
        }
        inventory = {
            level: _inventory(state, config, level)["capsule_calcium_lactate_mol"]
            for level, state in states.items()
        }
        self.assertAlmostEqual(inventory["0d"], inventory["1d"], places=12)
        self.assertAlmostEqual(inventory["0d"], inventory["2d"], places=12)

    def test_capsule_inventory_is_grid_independent(self):
        coarse = self.short_config()
        fine = copy.deepcopy(coarse)
        coarse.transport.nx_1d = 21
        fine.transport.nx_1d = 81
        _, coarse_profile = _capsule_profile_1d(coarse)
        _, fine_profile = _capsule_profile_1d(fine)
        coarse_inventory = _inventory(
            _initial_state(coarse, coarse_profile), coarse, "1d"
        )["capsule_calcium_lactate_mol"]
        fine_inventory = _inventory(
            _initial_state(fine, fine_profile), fine, "1d"
        )["capsule_calcium_lactate_mol"]
        self.assertAlmostEqual(coarse_inventory, fine_inventory, places=12)

    def test_spatial_calcite_summary_is_volume_integrated(self):
        config = self.short_config()
        config.transport.nx_1d = 9
        result = simulate_1d(config)
        final = result.frame.loc[result.frame["time_d"] == result.frame["time_d"].max()]
        integrated = float(np.sum(final["calcite_mol_m3"] * final["cell_volume_m3"]))
        self.assertAlmostEqual(
            integrated, result.diagnostics["final_inventory"]["calcite_mol"], places=12
        )

    def test_same_solid_volume_closes_narrower_crack_more(self):
        narrow = self.short_config()
        wide = copy.deepcopy(narrow)
        wide.transport.crack_width_mm = 2.0 * narrow.transport.crack_width_mm
        state = np.zeros((1, len(STATE_NAMES)))
        state[0, S["calcite_mol_m3"]] = 100.0
        volume = np.array([1.0e-9])
        wall_area = np.array([2.0e-6])
        narrow_metrics = repair_metrics(state, narrow, volume, wall_area)
        wide_metrics = repair_metrics(state, wide, volume, wall_area)
        self.assertGreater(
            narrow_metrics["crack_closure_ratio"][0],
            wide_metrics["crack_closure_ratio"][0],
        )

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

    def test_no_oxygen_prevents_calcite_mineralization(self):
        config = self.short_config()
        config.environment.oxygen_initial_mol_m3 = 0.0
        config.environment.oxygen_boundary_mol_m3 = 0.0
        config.environment.oxygen_transfer_s = 0.0
        result = simulate_0d(config)
        self.assertLess(result.summary["calcite_mol_m3_mean"], 1e-12)

    def test_response_delay_and_basal_leak_are_observable(self):
        delayed = self.short_config()
        delayed.simulation.days = 0.2
        delayed.kinetics.response_delay_h = 24.0
        immediate = copy.deepcopy(delayed)
        immediate.kinetics.response_delay_h = 0.0
        delayed_result = simulate_0d(delayed)
        immediate_result = simulate_0d(immediate)
        self.assertLessEqual(
            delayed_result.frame["activation_state"].iloc[-1],
            immediate_result.frame["activation_state"].iloc[-1] + 1e-9,
        )
        self.assertGreaterEqual(delayed_result.frame["false_activation_index"].iloc[0], 0.0)
        self.assertIn("premature_consumption_mol_m3", delayed_result.frame)
        self.assertIn("activation_delay_h", delayed_result.frame)

    def test_and_gate_requires_each_dynamic_signal(self):
        config = self.short_config()
        config.simulation.days = 0.5
        config.kinetics.gate_logic = "AND"
        config.kinetics.oxygen_rise_threshold_mol_m3_h = 1.0e6
        blocked_oxygen = simulate_0d(config)
        self.assertLess(blocked_oxygen.frame["environment_signal"].max(), 0.5)

        config = self.short_config()
        config.simulation.days = 0.2
        config.kinetics.gate_logic = "AND"
        config.kinetics.ph_drop_threshold_h = 1.0e6
        blocked_ph = simulate_0d(config)
        self.assertLess(blocked_ph.frame["environment_signal"].max(), 0.5)

    def test_and_gate_has_less_false_activation_than_or_gate(self):
        and_config = self.short_config()
        and_config.simulation.days = 0.2
        and_config.kinetics.gate_logic = "AND"
        or_config = copy.deepcopy(and_config)
        or_config.kinetics.gate_logic = "OR"
        and_result = simulate_0d(and_config)
        or_result = simulate_0d(or_config)
        self.assertLessEqual(
            and_result.frame["false_activation_index"].iloc[-1],
            or_result.frame["false_activation_index"].iloc[-1] + 1e-9,
        )

    def test_dry_condition_does_not_fully_activate_and_gate(self):
        config = self.short_config()
        config.environment.exposure = "dry"
        config.kinetics.gate_logic = "AND"
        result = simulate_0d(config)
        self.assertLess(result.frame["activation_state"].max(), 0.5)

    def test_and_gate_activates_after_duration_and_response_delay(self):
        config = self.short_config()
        config.simulation.days = 0.2
        config.environment.exposure = "continuous"
        config.kinetics.gate_logic = "AND"
        config.kinetics.oxygen_rise_threshold_mol_m3_h = -1.0
        config.kinetics.ph_drop_threshold_h = -1.0
        config.kinetics.ph_width = 100.0
        config.kinetics.oxygen_threshold_mol_m3 = 0.0
        config.kinetics.activation_duration_h = 0.1
        config.kinetics.response_delay_h = 0.1
        result = simulate_0d(config)
        self.assertGreater(result.frame["activation_state"].max(), 0.5)
        self.assertGreater(result.frame["activation_delay_h"].iloc[-1], 0.0)

    def test_basal_leak_increases_premature_consumption(self):
        low = self.short_config()
        low.simulation.days = 0.2
        low.kinetics.gate_logic = "AND"
        low.kinetics.basal_leak_fraction = 0.0
        low.kinetics.oxygen_rise_threshold_mol_m3_h = 1.0e6
        high = copy.deepcopy(low)
        high.kinetics.basal_leak_fraction = 0.1
        low_result = simulate_0d(low)
        high_result = simulate_0d(high)
        self.assertGreater(
            high_result.frame["premature_consumption_mol_m3"].iloc[-1],
            low_result.frame["premature_consumption_mol_m3"].iloc[-1],
        )

    def test_gate_logic_configuration_is_validated(self):
        config = self.short_config()
        config.kinetics.gate_logic = "XOR"
        with self.assertRaises(ValueError):
            config.validate()

    def test_precipitation_off_has_zero_closure_without_csh_fill(self):
        config = self.short_config()
        config.chemistry.calcite_rate_mol_m3_s = 0.0
        config.kinetics.capsule_csh_volume_fraction = 0.0
        result = simulate_0d(config)
        self.assertEqual(result.summary["mean_crack_closure_ratio"], 0.0)

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

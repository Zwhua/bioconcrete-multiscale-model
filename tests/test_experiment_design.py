import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from bioconcrete.experiment_design import (
    d_optimal_score,
    fisher_information,
    greedy_d_optimal,
    numerical_candidate_table,
    rank_experiments,
)


def analytic_trajectory(config):
    times = np.array([0.0, 1.0, 7.0, 14.0, 21.0, 28.0])
    activity = config.kinetics.activity_multiplier
    release = config.kinetics.capsule_release_s / 1.3e-6
    calcite = config.chemistry.calcite_rate_mol_m3_s / 1e-4
    wall = config.chemistry.wall_deposition_fraction
    csh = config.kinetics.capsule_csh_volume_fraction
    return pd.DataFrame({
        "time_d": times,
        "crack_closure_ratio": times / 28.0 * (activity + wall + csh),
        "calcite_mol_m3": times * calcite,
        "lactate_mol_m3": times * release,
        "oxygen_mol_m3": np.full(len(times), config.environment.oxygen_transfer_s * 1e5),
        "calcium_mol_m3": times * config.chemistry.portlandite_mol_m3 / 5000.0,
        "ph": np.full(len(times), config.environment.ph),
        "csh_volume_fraction": times / 28.0 * csh,
    })


class ExperimentDesignTests(unittest.TestCase):
    def test_independent_measurements_have_more_information(self):
        duplicate = np.array([[1, 0], [1, 0]], dtype=float)
        independent = np.eye(2)
        self.assertGreater(d_optimal_score(independent), d_optimal_score(duplicate))

    def test_measurement_noise_reduces_information(self):
        jacobian = np.array([[1.0, 2.0]])
        quiet = fisher_information(jacobian, np.array([[0.1 ** 2]]))
        noisy = fisher_information(jacobian, np.array([[1.0]]))
        self.assertGreater(np.trace(quiet), np.trace(noisy))

    def test_numerical_sensitivity_comes_from_runner(self):
        frame = numerical_candidate_table(
            trajectory_runner=analytic_trajectory, widths=(0.3,), doses=(1.0,),
            wettings=(12.0,), times=(7.0,), systems=("complete",),
        )
        activity = frame.loc[frame["observable"] == "crack_closure_ratio",
                             "sensitivity_kinetics.activity_multiplier"].iloc[0]
        self.assertGreater(activity, 0.0)

    def test_greedy_selection_is_complementary(self):
        rows = []
        for identifier, vector in (("a", [1, 0, 0, 0, 0]), ("b", [0, 1, 0, 0, 0]),
                                   ("duplicate", [1, 0, 0, 0, 0])):
            row = {"experiment_id": identifier, "measurement_sd": 1.0, "execution_cost": 1.0}
            row.update({"sensitivity_{}".format(name): value for name, value in zip(
                ("kinetics.capsule_release_s", "kinetics.activity_multiplier",
                 "chemistry.calcite_rate_mol_m3_s", "chemistry.wall_deposition_fraction",
                 "kinetics.capsule_csh_volume_fraction"), vector)})
            rows.append(row)
        selected, _ = greedy_d_optimal(pd.DataFrame(rows), 2)
        self.assertEqual(set(selected["experiment_id"]), {"a", "b"})

    def test_ranked_plan_has_required_top_experiments(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = rank_experiments(Path(directory), trajectory_runner=analytic_trajectory)
            frame = pd.read_csv(Path(directory) / "recommended_experiments.csv")
            self.assertEqual(summary["minimum_executable_count"], 5)
            self.assertEqual(len(frame), 10)
            self.assertFalse(summary["experimental_validation"])
            self.assertEqual(summary["method"], "numerical greedy D-optimal")


if __name__ == "__main__":
    unittest.main()

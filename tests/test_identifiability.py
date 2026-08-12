import tempfile
import unittest
from pathlib import Path

import numpy as np

from bioconcrete.config import ModelConfig
from bioconcrete.identifiability import identifiability_analysis, matrix_diagnostics


class IdentifiabilityTests(unittest.TestCase):
    def test_matrix_diagnostics_detects_collinear_parameters(self):
        matrix = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
        table, correlation, summary = matrix_diagnostics(matrix, ["a", "b"])
        self.assertEqual(summary["fim_rank"], 1)
        self.assertEqual(summary["unidentifiable_combinations"], 1)
        self.assertTrue((table["identifiability_label"] != "estimable").all())
        self.assertAlmostEqual(abs(correlation.loc["a", "b"]), 1.0)

    def test_matrix_diagnostics_accepts_independent_parameters(self):
        table, _, summary = matrix_diagnostics(np.eye(3), ["a", "b", "c"])
        self.assertTrue(summary["full_rank"])
        self.assertTrue((table["identifiability_label"] == "estimable").all())

    def test_workflow_writes_required_outputs(self):
        config = ModelConfig()
        config.simulation.days = 0.1
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary = identifiability_analysis(
                output, config,
                parameters=("kinetics.activity_multiplier",),
                outputs=("crack_closure_ratio",), times_d=(0.05, 0.1),
            )
            self.assertFalse(summary["calibrated"])
            self.assertEqual(summary["profile_likelihood"], "not_executable_without_observed_likelihood")
            for name in (
                "local_sensitivity.csv", "parameter_correlation.csv",
                "identifiability_table.csv", "recommended_measurements.csv",
                "identifiability_summary.json",
            ):
                self.assertTrue((output / name).exists())


if __name__ == "__main__":
    unittest.main()

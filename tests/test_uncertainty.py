import tempfile
import unittest
from pathlib import Path

import pandas as pd

from bioconcrete.config import ModelConfig
from bioconcrete.uncertainty import prior_predictive, sample_prior


class UncertaintyTests(unittest.TestCase):
    def test_uniform_and_log_uniform_priors(self):
        self.assertAlmostEqual(sample_prior({"distribution": "uniform", "low": 2, "high": 4}, 0.5), 3.0)
        self.assertAlmostEqual(sample_prior({"distribution": "log_uniform", "low": 1, "high": 100}, 0.5), 10.0)

    def test_prior_predictive_is_reproducible_and_resumable(self):
        config = ModelConfig()
        config.simulation.days = 0.05
        priors = {"kinetics.activity_multiplier": {"distribution": "uniform", "low": 0.8, "high": 1.2}}
        scenarios = {"transport.crack_width_mm": {"distribution": "uniform", "low": 0.2, "high": 0.4}}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            first = prior_predictive(output, config, samples=2, seed=7, priors=priors, scenarios=scenarios)
            original = pd.read_csv(output / "prior_predictive_samples.csv")
            second = prior_predictive(output, config, samples=2, seed=7, priors=priors, scenarios=scenarios, resume=True)
            resumed = pd.read_csv(output / "prior_predictive_samples.csv")
            pd.testing.assert_frame_equal(original, resumed)
            self.assertEqual(first["interval_type"], "prior_predictive_interval")
            self.assertFalse(second["calibrated_prediction_interval"])
            self.assertEqual(set(original["evidence_class"]), {"prior_predictive_model_output"})


if __name__ == "__main__":
    unittest.main()

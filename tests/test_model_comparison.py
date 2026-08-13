import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from bioconcrete.config import ModelConfig
from bioconcrete.model_comparison import compare_structures, information_criteria


class ModelComparisonTests(unittest.TestCase):
    def test_information_criteria_formula(self):
        observed = np.array([0, 1, 2, 3, 4, 5], float)
        predicted = observed + 0.1
        result = information_criteria(observed, predicted, 2)
        rss = np.sum((predicted - observed) ** 2)
        self.assertAlmostEqual(result["aic"], 6 * np.log(rss / 6) + 4)
        self.assertAlmostEqual(result["mae"], 0.1)

    def test_no_data_means_no_aic(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = compare_structures(Path(directory), ModelConfig())
            frame = pd.read_csv(Path(directory) / "model_structure_comparison.csv")
            self.assertFalse(summary["information_criteria_available"])
            self.assertTrue(frame["aic"].isna().all())

    def test_predictions_use_actual_times_and_conditions(self):
        observations = pd.DataFrame({
            "specimen_id": ["A", "A", "B", "B", "C", "D"],
            "split": ["train", "train", "train", "train", "internal_test", "internal_test"],
            "time_d": [1.0, 7.0, 1.0, 7.0, 14.0, 28.0],
            "initial_crack_width_mm": [0.1, 0.1, 0.5, 0.5, 0.3, 0.3],
            "crack_closure_ratio": [0.01, 0.07, 0.005, 0.03, 0.12, 0.20],
        })
        calls = []

        def predictor(config, frame):
            calls.append(frame[["time_d", "initial_crack_width_mm"]].copy())
            return pd.DataFrame({
                "crack_closure_ratio": frame["time_d"].to_numpy(float) / 100.0
                / frame["initial_crack_width_mm"].to_numpy(float)
            })

        with tempfile.TemporaryDirectory() as directory:
            summary = compare_structures(Path(directory), observations=observations, predictor=predictor)
            predictions = pd.read_csv(Path(directory) / "model_structure_predictions.csv")
            full = predictions.loc[predictions["structure"] == "full_mechanistic"]
            self.assertEqual(full["time_d"].tolist(), observations["time_d"].tolist())
            self.assertGreater(full["predicted"].nunique(), 1)
            self.assertEqual(len(calls), 4)
            self.assertTrue(summary["condition_resolved"])

    def test_unreviewed_rows_do_not_gain_calibrated_label(self):
        observations = pd.DataFrame({
            "specimen_id": ["A", "A", "B", "B"], "time_d": [1, 7, 1, 7],
            "crack_closure_ratio": [0.01, 0.1, 0.02, 0.12],
            "curation_status": ["candidate_only"] * 4,
        })
        predictor = lambda config, frame: pd.DataFrame({"crack_closure_ratio": np.zeros(len(frame))})
        with tempfile.TemporaryDirectory() as directory:
            compare_structures(Path(directory), observations=observations, predictor=predictor)
            metrics = pd.read_csv(Path(directory) / "model_structure_comparison.csv")
            self.assertEqual(set(metrics["evidence_label"]), {"uncalibrated prospective model prediction"})


if __name__ == "__main__":
    unittest.main()

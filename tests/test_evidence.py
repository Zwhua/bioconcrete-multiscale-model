import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from bioconcrete.config import ModelConfig
from bioconcrete.evidence import (
    _digest_config,
    _fit_empirical_closure,
    _metrics,
    _sigma_for_output,
    _standardized_residuals,
    calibrate_public,
)
from bioconcrete.public_data import OBSERVATION_COLUMNS, grouped_split, prepare_public_data


class EvidenceTests(unittest.TestCase):
    def test_grouped_split_never_splits_one_specimen(self):
        specimens = ["A", "A", "B", "B", "C"]
        first = grouped_split(specimens)
        second = grouped_split(reversed(specimens))
        self.assertEqual(first, second)
        self.assertEqual(first["A"], first["A"])

    def test_normalized_schema_is_stable_for_empty_download(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "raw" / "empty").mkdir(parents=True)
            summary = prepare_public_data("empty", root)
            frame = pd.read_csv(root / "derived" / "empty" / "observations.csv")
            self.assertEqual(list(frame.columns), OBSERVATION_COLUMNS)
            self.assertEqual(summary["rows"], 0)
            self.assertFalse(summary["calibration_ready"])

    def test_frozen_configuration_digest_detects_changes(self):
        config = ModelConfig()
        original = _digest_config(config)
        config.chemistry.calcite_rate_mol_m3_s *= 2.0
        self.assertNotEqual(original, _digest_config(config))

    def test_metrics_use_correct_aic_and_aicc(self):
        observed = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        predicted = np.array([0.1, 0.9, 2.2, 2.8, 4.1, 4.9])
        metrics = _metrics(observed, predicted, parameter_count=2)
        rss = float(np.sum((predicted - observed) ** 2))
        expected_aic = len(observed) * np.log(rss / len(observed)) + 4.0
        expected_aicc = expected_aic + 2.0 * 2 * 3 / (len(observed) - 2 - 1)
        self.assertAlmostEqual(metrics["aic"], expected_aic)
        self.assertAlmostEqual(metrics["aicc"], expected_aicc)
        self.assertIn("calibration_slope", metrics)
        self.assertIn("calibration_intercept", metrics)

    def test_sigma_priority_uses_provided_measurement_sd(self):
        frame = pd.DataFrame({
            "initial_crack_width_mm": [0.2, 0.2, 0.2],
            "measurement_sd": [0.002, np.nan, np.nan],
        })
        observed = np.array([0.1, 0.2, 0.3])
        sigma, audit = _sigma_for_output(
            frame, "crack_closure_ratio", observed,
            independent_error={"crack_closure_ratio": 0.03},
        )
        self.assertAlmostEqual(sigma[0], 0.01)
        np.testing.assert_allclose(sigma[1:], 0.03)
        self.assertIn("provided:measurement_sd_propagated", audit["sources"])
        self.assertIn("independent_measurement_error", audit["sources"])

    def test_standardized_residuals_do_not_mix_physical_units(self):
        frame = pd.DataFrame({"time_d": [1.0, 2.0], "specimen_id": ["A", "B"]})
        predictions = pd.DataFrame({
            "lactate_mM": [12.0, 18.0],
            "calcite_mass_mg": [1.5, 2.5],
            "crack_closure_ratio": [0.0, 0.0],
            "permeability_ratio": [0.0, 0.0],
            "ph": [0.0, 0.0], "activation_state": [0.0, 0.0],
            "cumulative_activity_h": [0.0, 0.0],
        })
        observed_frame = frame.assign(lactate_mM=[10.0, 20.0], calcite_mass_mg=[1.0, 3.0])
        sigma = {"lactate_mM": np.array([2.0, 2.0]), "calcite_mass_mg": np.array([0.5, 0.5])}
        with patch("bioconcrete.evidence._predict_outputs", return_value=predictions):
            residuals, contributions, _ = _standardized_residuals(
                ModelConfig(), observed_frame, ["lactate_mM", "calcite_mass_mg"], sigma
            )
        np.testing.assert_allclose(residuals, [1.0, -1.0, 1.0, -1.0])
        self.assertAlmostEqual(contributions["lactate_mM"], 2.0)
        self.assertAlmostEqual(contributions["calcite_mass_mg"], 2.0)

    def test_empirical_baseline_is_fitted_from_training_closure(self):
        times = np.array([0.0, 2.0, 5.0, 10.0, 20.0, 28.0])
        expected = np.array([0.0, 0.04, 0.09, 0.15, 0.22, 0.25])
        frame = pd.DataFrame({"time_d": times, "crack_closure_ratio": expected})
        parameters, artifact = _fit_empirical_closure(frame, np.full(len(frame), 0.01))
        self.assertIsNotNone(parameters)
        self.assertEqual(artifact["status"], "fitted_on_training_closure")
        self.assertGreater(artifact["h_inf"], 0.0)
        self.assertGreater(artifact["rate_d_inv"], 0.0)

    def test_no_calcite_observation_keeps_precipitation_rate_at_prior(self):
        frame = pd.DataFrame({
            "dataset_id": ["public_A"] * 8,
            "specimen_id": ["A"] * 4 + ["B"] * 4,
            "split": ["train"] * 8,
            "time_d": [0.0, 1.0, 2.0, 3.0] * 2,
            "lactate_mM": [0.0, 1.0, 2.0, 3.0] * 2,
            "crack_closure_ratio": [0.0, 0.01, 0.02, 0.03] * 2,
            "permeability_ratio": [1.0, 0.99, 0.98, 0.97] * 2,
        })

        def fake_predictions(config, observations):
            time = pd.to_numeric(observations["time_d"]).to_numpy(float)
            release = config.kinetics.capsule_release_s / 1.3e-6
            activity = config.kinetics.activity_multiplier
            wall = config.chemistry.wall_deposition_fraction
            return pd.DataFrame({
                "lactate_mM": release * activity * time,
                "calcite_mass_mg": np.zeros(len(time)),
                "crack_closure_ratio": wall * time / 75.0,
                "permeability_ratio": 1.0 - wall * time / 75.0,
                "ph": np.full(len(time), 11.5),
                "activation_state": np.full(len(time), activity / 5.0),
                "cumulative_activity_h": activity * time,
            })

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data_path = root / "observations.csv"
            frame.to_csv(data_path, index=False)
            prior = ModelConfig().chemistry.calcite_rate_mol_m3_s
            with patch("bioconcrete.evidence._predict_outputs", side_effect=fake_predictions):
                result = calibrate_public(data_path, root / "run", bootstrap_samples=2)
            fitted = ModelConfig.load(root / "run" / "frozen_config.json")
            self.assertAlmostEqual(fitted.chemistry.calcite_rate_mol_m3_s, prior)
            self.assertIn("chemistry.calcite_rate_mol_m3_s", result["fixed_to_prior"])
            self.assertTrue((root / "run" / "calibration_metrics.json").exists())
            self.assertTrue((root / "run" / "train_predictions.csv").exists())
            self.assertTrue((root / "run" / "calibrated_prediction_intervals.csv").exists())
            self.assertTrue((root / "run" / "calibration_plots" / "lactate_mM_observed_vs_predicted.png").exists())
            intervals = pd.read_csv(root / "run" / "calibrated_prediction_intervals.csv")
            self.assertEqual(set(intervals["sample_count"]), {2})
            predictions = pd.read_csv(root / "run" / "train_predictions.csv")
            self.assertTrue(predictions["prediction_low_95"].notna().all())
            metrics = pd.read_json(root / "run" / "calibration_metrics.json", typ="series")
            self.assertIn("per_output", metrics)

    def test_closure_only_data_cannot_fit_upstream_parameters(self):
        frame = pd.DataFrame({
            "dataset_id": ["public_A"] * 6,
            "specimen_id": ["A", "A", "B", "B", "C", "C"],
            "split": ["train"] * 6,
            "time_d": [1.0, 7.0] * 3,
            "crack_closure_ratio": [0.01, 0.12, 0.02, 0.14, 0.01, 0.13],
        })

        def fake_predictions(config, observations):
            time = pd.to_numeric(observations["time_d"]).to_numpy(float)
            release = config.kinetics.capsule_release_s / 1.3e-6
            activity = config.kinetics.activity_multiplier
            calcite = config.chemistry.calcite_rate_mol_m3_s / 2.5e-5
            wall = config.chemistry.wall_deposition_fraction
            closure = release * activity * calcite * wall * time / 120.0
            return pd.DataFrame({
                "lactate_mM": np.zeros(len(time)),
                "calcite_mass_mg": np.zeros(len(time)),
                "crack_closure_ratio": closure,
                "permeability_ratio": 1.0 - closure,
                "ph": np.full(len(time), 11.5),
                "activation_state": np.zeros(len(time)),
                "cumulative_activity_h": np.zeros(len(time)),
            })

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data_path = root / "observations.csv"
            frame.to_csv(data_path, index=False)
            prior = ModelConfig()
            with patch("bioconcrete.evidence._predict_outputs", side_effect=fake_predictions):
                result = calibrate_public(data_path, root / "run", bootstrap_samples=0)
            fitted = ModelConfig.load(root / "run" / "frozen_config.json")
            self.assertAlmostEqual(fitted.kinetics.capsule_release_s, prior.kinetics.capsule_release_s)
            self.assertAlmostEqual(fitted.kinetics.activity_multiplier, prior.kinetics.activity_multiplier)
            self.assertAlmostEqual(
                fitted.chemistry.calcite_rate_mol_m3_s,
                prior.chemistry.calcite_rate_mol_m3_s,
            )
            self.assertEqual(result["fitted_parameters"], ["chemistry.wall_deposition_fraction"])

    def test_candidate_extraction_is_rejected_as_calibration_evidence(self):
        frame = pd.DataFrame({
            "dataset_id": ["public_A"], "specimen_id": ["A"],
            "split": ["train"], "time_d": [7.0],
            "crack_closure_ratio": [0.2], "curation_status": ["candidate_only"],
        })
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data_path = root / "candidate.csv"
            frame.to_csv(data_path, index=False)
            with self.assertRaisesRegex(ValueError, "Candidate extraction"):
                calibrate_public(data_path, root / "run", bootstrap_samples=0)


if __name__ == "__main__":
    unittest.main()

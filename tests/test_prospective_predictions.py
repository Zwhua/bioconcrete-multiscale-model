import unittest
from pathlib import Path

import yaml


class ProspectivePredictionTests(unittest.TestCase):
    def test_predictions_are_falsifiable_and_not_claimed_validated(self):
        root = Path(__file__).resolve().parents[1]
        content = yaml.safe_load((root / "PROSPECTIVE_PREDICTIONS.yml").read_text(encoding="utf-8"))
        self.assertFalse(content["team_wet_lab_data_available"])
        self.assertGreaterEqual(len(content["predictions"]), 3)
        required = {
            "prediction_id", "input_conditions", "point_prediction",
            "success_threshold", "failure_threshold", "required_future_measurement",
            "evidence_class", "falsifies",
        }
        for prediction in content["predictions"]:
            self.assertTrue(required.issubset(prediction))
            self.assertIn("requiring future experimental validation", prediction["evidence_class"])


if __name__ == "__main__":
    unittest.main()

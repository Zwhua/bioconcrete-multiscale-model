import tempfile
import unittest
from pathlib import Path

import pandas as pd

from bioconcrete.decision_support import REQUIRED_FIELDS, generate_decision_support


class DecisionSupportTests(unittest.TestCase):
    def test_decision_tables_contain_evidence_and_thresholds(self):
        frame = pd.DataFrame([{
            "crack_width_mm": 0.3, "agent_dosage": 1.0, "wet_hours_per_day": 12,
            "activity_multiplier": 1.0, "response_delay_h": 4,
            "basal_leak_fraction": 0.01, "closure_28d": 0.2,
            "decision": "recommended",
        }, {
            "crack_width_mm": 0.5, "agent_dosage": 0.5, "wet_hours_per_day": 6,
            "activity_multiplier": 0.5, "response_delay_h": 24,
            "basal_leak_fraction": 0.1, "closure_28d": 0.01,
            "decision": "not_recommended",
        }])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "design.csv"
            frame.to_csv(source, index=False)
            generate_decision_support(source, root / "out", "config", "code")
            recommended = pd.read_csv(root / "out" / "recommended_conditions.csv")
            rejected = pd.read_csv(root / "out" / "rejected_conditions.csv")
            self.assertTrue(set(REQUIRED_FIELDS).issubset(recommended.columns))
            self.assertEqual(len(rejected), 1)
            self.assertEqual(set(recommended["evidence_level"]), {"uncalibrated prospective model prediction"})


if __name__ == "__main__":
    unittest.main()

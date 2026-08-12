import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from bioconcrete.design import _evaluate_scenario, design_matrix


def fake_scenario(payload):
    setting, _ = payload
    value = sum(float(item) for item in setting.values())
    return {**setting, "scenario_id": __import__("bioconcrete.design", fromlist=["_scenario_id"])._scenario_id(setting),
            "closure_28d": value / 100, "time_to_50pct_d": float("nan"),
            "permeability_ratio": 1 - value / 1000, "closure_per_agent": value / 100,
            "premature_consumption": 0.1, "target_probability": float("nan"),
            "runtime_s": 0.0, "dominant_bottleneck": "geometry_limited", "bottleneck_score": 0.5}


class DesignParallelTests(unittest.TestCase):
    def test_resume_does_not_duplicate_completed_scenarios(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch("bioconcrete.design._evaluate_scenario", side_effect=fake_scenario):
                design_matrix(root / "PREREGISTERED_SCENARIOS.yml", output, limit=2)
                first = pd.read_csv(output / "design_matrix.csv")
                design_matrix(root / "PREREGISTERED_SCENARIOS.yml", output, limit=2, resume=True)
                second = pd.read_csv(output / "design_matrix.csv")
            self.assertEqual(len(first), 2)
            self.assertEqual(len(second), 2)
            self.assertEqual(second["scenario_id"].nunique(), 2)


if __name__ == "__main__":
    unittest.main()

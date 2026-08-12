import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from bioconcrete.experiment_design import d_optimal_score, rank_experiments


class ExperimentDesignTests(unittest.TestCase):
    def test_independent_measurements_have_more_information(self):
        duplicate = np.array([[1, 0], [1, 0]], dtype=float)
        independent = np.eye(2)
        self.assertGreater(d_optimal_score(independent), d_optimal_score(duplicate))

    def test_ranked_plan_has_required_top_experiments(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = rank_experiments(Path(directory))
            frame = pd.read_csv(Path(directory) / "recommended_experiments.csv")
            self.assertEqual(summary["minimum_executable_count"], 5)
            self.assertEqual(len(frame), 10)
            self.assertFalse(summary["experimental_validation"])


if __name__ == "__main__":
    unittest.main()

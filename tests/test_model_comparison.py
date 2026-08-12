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

    def test_no_data_means_no_aic(self):
        config = ModelConfig(); config.simulation.days = 0.05
        with tempfile.TemporaryDirectory() as directory:
            summary = compare_structures(Path(directory), config)
            frame = pd.read_csv(Path(directory) / "model_structure_comparison.csv")
            self.assertFalse(summary["information_criteria_available"])
            self.assertTrue(frame["aic"].isna().all())


if __name__ == "__main__":
    unittest.main()

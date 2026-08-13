import tempfile
import unittest
from pathlib import Path

import pandas as pd

from bioconcrete.biological_design import generate_biological_design


class BiologicalDesignTests(unittest.TestCase):
    def test_mapping_is_anonymous_and_unit_complete(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            summary = generate_biological_design(output)
            mapping = pd.read_csv(output / "part_parameter_mapping.csv")
            self.assertFalse(summary["contains_sequences"])
            self.assertTrue(mapping["unit"].notna().all())
            self.assertTrue(mapping["required_measurement"].notna().all())
            forbidden = "sequence mutation plasmid vector construction"
            text = " ".join(mapping.astype(str).to_numpy().ravel()).lower()
            for word in forbidden.split():
                self.assertNotIn(word, text)

    def test_candidate_scores_are_not_measured_activity(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            generate_biological_design(output)
            predictions = pd.read_csv(output / "construct_predictions.csv")
            self.assertTrue(predictions["prediction_type"].str.contains("not measured").all())


if __name__ == "__main__":
    unittest.main()

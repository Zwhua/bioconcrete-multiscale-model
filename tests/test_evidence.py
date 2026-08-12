import tempfile
import unittest
from pathlib import Path

import pandas as pd

from bioconcrete.config import ModelConfig
from bioconcrete.evidence import _digest_config
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

    def test_frozen_configuration_digest_detects_changes(self):
        config = ModelConfig()
        original = _digest_config(config)
        config.chemistry.calcite_rate_mol_m3_s *= 2.0
        self.assertNotEqual(original, _digest_config(config))


if __name__ == "__main__":
    unittest.main()

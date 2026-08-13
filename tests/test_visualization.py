import json
import tempfile
import unittest
from pathlib import Path

from bioconcrete.visualization import render_figures


class VisualizationTests(unittest.TestCase):
    def test_missing_inputs_are_reported_without_substitution(self):
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory) / "run"; run.mkdir()
            result = render_figures(run)
            self.assertIn("figure01_model_to_decision", result["generated"])
            self.assertIn("figure08_multiscale", result["missing"])
            self.assertFalse(result["substitute_data_generated"])
            self.assertTrue((run / "figures" / "figure01_model_to_decision.png").exists())
            self.assertTrue((run / "figures" / "figure01_model_to_decision.svg").exists())
            manifest = json.loads((run / "figures" / "figure_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["missing"], result["missing"])


if __name__ == "__main__":
    unittest.main()

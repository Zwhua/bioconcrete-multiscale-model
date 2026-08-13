import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioconcrete.release_analysis import RELEASE_DIRECTORIES, release_analysis


class ReleaseArtifactTests(unittest.TestCase):
    def test_initialize_never_claims_formal_completion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "PREREGISTERED_SCENARIOS.yml").write_text("factors: {}", encoding="utf-8")
            with patch("bioconcrete.manifest._git", return_value="clean-test"):
                result = release_analysis(root, initialize_only=True)
            release = root / "model_runs" / "v0.5.0"
            self.assertFalse(result["formal_analyses_complete"])
            for name in RELEASE_DIRECTORIES:
                self.assertTrue((release / name).is_dir())
            manifest = json.loads((release / "release_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["evidence_label"], "uncalibrated prospective model prediction")
            self.assertEqual(manifest["team_wet_lab_rows"], 0)


if __name__ == "__main__":
    unittest.main()

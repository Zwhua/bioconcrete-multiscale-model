import tempfile
import unittest
from pathlib import Path

from bioconcrete.dashboard import generate_dashboard


class DashboardTests(unittest.TestCase):
    def test_dashboard_marks_missing_evidence(self):
        root = Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as directory:
            result = generate_dashboard(root, Path(directory), Path(directory) / "missing-run")
            content = Path(result["dashboard"]).read_text(encoding="utf-8")
            self.assertIn("Team wet-lab observations: 0", content)
            self.assertFalse(result["generated_data"])
            self.assertIn("What did the model change?", content)
            self.assertTrue(Path(result["static_fallback"]).exists())


if __name__ == "__main__":
    unittest.main()

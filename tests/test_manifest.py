import tempfile
import unittest
from pathlib import Path

from bioconcrete.config import ModelConfig
from bioconcrete.manifest import create_manifest, finish_manifest, sha256_file


class ManifestTests(unittest.TestCase):
    def test_manifest_contains_required_provenance(self):
        root = Path(__file__).resolve().parents[1]
        manifest = create_manifest(root, ModelConfig(), ["test"], 2026)
        required = {
            "model_version", "git_commit", "git_worktree_dirty", "config_sha256",
            "input_sha256", "random_seed", "started_at_utc", "completed_at_utc",
            "software_environment", "command", "status",
        }
        self.assertTrue(required.issubset(manifest))
        completed = finish_manifest(manifest)
        self.assertEqual(completed["status"], "complete")
        self.assertIsNotNone(completed["completed_at_utc"])

    def test_file_digest_is_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.txt"
            path.write_text("evidence", encoding="utf-8")
            self.assertEqual(sha256_file(path), sha256_file(path))


if __name__ == "__main__":
    unittest.main()

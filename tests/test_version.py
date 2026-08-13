import re
import unittest
from pathlib import Path

import yaml

import bioconcrete


class VersionTests(unittest.TestCase):
    def test_release_version_is_consistent(self):
        root = Path(__file__).resolve().parents[1]
        pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
        project_version = re.search(
            r'^version = "([^"]+)"$', pyproject, re.MULTILINE
        ).group(1)
        citation = yaml.safe_load((root / "CITATION.cff").read_text(encoding="utf-8"))
        readme = (root / "README.md").read_text(encoding="utf-8")
        readme_zh = (root / "README.zh-CN.md").read_text(encoding="utf-8")
        self.assertEqual(bioconcrete.__version__, project_version)
        self.assertEqual(str(citation["version"]), project_version)
        self.assertIn("version-v{}".format(project_version), readme)
        self.assertIn("Current release: [v{}]".format(project_version), readme)
        self.assertIn("当前版本：[v{}]".format(project_version), readme_zh)


if __name__ == "__main__":
    unittest.main()

"""Repository-homepage consistency checks."""

from pathlib import Path
import re
import shlex
import unittest

from bioconcrete.cli import build_parser


ROOT = Path(__file__).resolve().parents[1]
READMES = (ROOT / "README.md", ROOT / "README.zh-CN.md")


class ReadmeTests(unittest.TestCase):
    def test_local_markdown_targets_exist(self) -> None:
        pattern = re.compile(r"!?(?:\[[^]]*\])\(([^)]+)\)")
        for readme in READMES:
            for target in pattern.findall(readme.read_text(encoding="utf-8")):
                target = target.split("#", 1)[0]
                if not target or "://" in target or target.startswith("#"):
                    continue
                self.assertTrue((ROOT / target).exists(), "{} -> {}".format(readme.name, target))

    def test_bilingual_evidence_and_numbers_match(self) -> None:
        for readme in READMES:
            text = readme.read_text(encoding="utf-8")
            for required in ("v0.5.1", "2.0819%", "2.08023", "0.00168", "0.9388", "82"):
                self.assertIn(required, text)
            self.assertNotIn("public-data-supported prediction", text)
            self.assertNotIn("fully validated", text.lower())

    def test_current_readme_does_not_present_v4_commands(self) -> None:
        for readme in READMES:
            text = readme.read_text(encoding="utf-8")
            self.assertNotIn("V4 analysis commands", text)
            self.assertNotIn("V4 分析命令", text)

    def test_bioconcrete_commands_are_parseable(self) -> None:
        parser = build_parser()
        command = re.compile(r"^python -m bioconcrete (.+)$", re.MULTILINE)
        for readme in READMES:
            for arguments in command.findall(readme.read_text(encoding="utf-8")):
                parsed = parser.parse_args(shlex.split(arguments))
                self.assertIsNotNone(parsed.command)


if __name__ == "__main__":
    unittest.main()

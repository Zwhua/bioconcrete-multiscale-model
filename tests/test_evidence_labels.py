import tempfile
import unittest
from pathlib import Path

from bioconcrete.evidence_state import (
    EXTERNALLY_EVALUATED,
    PUBLIC_CALIBRATED,
    UNCALIBRATED,
    evidence_state,
    validate_evidence_label,
)


class EvidenceLabelTests(unittest.TestCase):
    def test_state_cannot_skip_missing_calibration(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run, external = root / "run", root / "external"
            run.mkdir()
            external.mkdir()
            self.assertEqual(evidence_state(run, external).label, UNCALIBRATED)
            (external / "external_validation.json").write_text("{}", encoding="utf-8")
            self.assertEqual(evidence_state(run, external).label, UNCALIBRATED)
            (run / "frozen_run.json").write_text("{}", encoding="utf-8")
            self.assertEqual(evidence_state(run, external).label, EXTERNALLY_EVALUATED)

    def test_calibration_state_requires_frozen_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "frozen_run.json").write_text("{}", encoding="utf-8")
            self.assertEqual(evidence_state(root, root / "missing").label, PUBLIC_CALIBRATED)

    def test_ad_hoc_supported_claim_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_evidence_label("public-data-supported prediction")


if __name__ == "__main__":
    unittest.main()

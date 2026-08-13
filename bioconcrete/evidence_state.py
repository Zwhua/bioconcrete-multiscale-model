"""Central evidence-state vocabulary and claim gating."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict


UNCALIBRATED = "uncalibrated prospective model prediction"
PUBLIC_CALIBRATED = "public-data-calibrated prospective prediction"
EXTERNALLY_EVALUATED = "externally evaluated"
TEAM_OBSERVATION = "team wet-lab observation"
ALLOWED_EVIDENCE_LABELS = (
    UNCALIBRATED, PUBLIC_CALIBRATED, EXTERNALLY_EVALUATED, TEAM_OBSERVATION,
)


@dataclass(frozen=True)
class EvidenceState:
    """Auditable state derived only from completed evidence artifacts."""

    label: str
    public_calibration_complete: bool
    external_evaluation_complete: bool
    team_wet_lab_rows: int = 0

    def to_dict(self) -> Dict[str, object]:
        return {
            "evidence_label": self.label,
            "public_calibration_complete": self.public_calibration_complete,
            "external_evaluation_complete": self.external_evaluation_complete,
            "team_wet_lab_rows": self.team_wet_lab_rows,
        }


def evidence_state(run_dir: Path, external_dir: Path, team_wet_lab_rows: int = 0) -> EvidenceState:
    """Resolve the strongest completed state without inferring missing evidence."""

    calibrated = (run_dir / "frozen_run.json").exists()
    externally_evaluated = calibrated and (external_dir / "external_validation.json").exists()
    if team_wet_lab_rows > 0:
        label = TEAM_OBSERVATION
    elif externally_evaluated:
        label = EXTERNALLY_EVALUATED
    elif calibrated:
        label = PUBLIC_CALIBRATED
    else:
        label = UNCALIBRATED
    return EvidenceState(label, calibrated, externally_evaluated, team_wet_lab_rows)


def validate_evidence_label(label: str) -> str:
    """Reject ad-hoc labels that can silently overstate model evidence."""

    if label not in ALLOWED_EVIDENCE_LABELS:
        raise ValueError("Unsupported evidence label: {}".format(label))
    return label

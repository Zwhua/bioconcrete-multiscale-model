"""V0.5.0 release-analysis orchestration and artifact completeness checks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

from .biological_design import generate_biological_design
from .config import ModelConfig
from .evidence_state import UNCALIBRATED
from .manifest import create_manifest, finish_manifest, write_manifest


RELEASE_DIRECTORIES = (
    "baseline", "validation", "identifiability", "uncertainty", "sensitivity",
    "design_matrix", "model_comparison", "counterfactual_bottleneck",
    "experiment_design", "biological_design", "decision_support", "dashboard",
)


def release_analysis(
    project_root: Path, version: str = "0.5.0", config: Optional[ModelConfig] = None,
    workers: int = 16, resume: bool = False, initialize_only: bool = False,
) -> Dict[str, object]:
    """Initialize a release run without claiming unexecuted long analyses."""

    if version != "0.5.0":
        raise ValueError("This orchestrator is pinned to the v0.5.0 evidence contract")
    base = project_root / "model_runs" / "v0.5.0"
    base.mkdir(parents=True, exist_ok=True)
    for name in RELEASE_DIRECTORIES:
        (base / name).mkdir(exist_ok=True)
    model_config = config or ModelConfig()
    manifest = create_manifest(
        project_root, model_config, ["release-analysis", "--version", version], 2026,
        {"preregister": project_root / "PREREGISTERED_SCENARIOS.yml"},
        status="initialized" if initialize_only else "incomplete",
    )
    manifest.update({
        "evidence_label": UNCALIBRATED, "workers": workers, "resume": resume,
        "team_wet_lab_rows": 0, "public_calibration_complete": False,
        "external_evaluation_complete": False,
        "required_directories": list(RELEASE_DIRECTORIES),
        "formal_result_policy": "missing analyses remain missing; smoke output is never promoted",
    })
    generate_biological_design(base / "biological_design")
    status = {
        "release_version": version, "release_directory": str(base),
        "initialized": True, "formal_analyses_complete": False,
        "long_run_required": True,
        "next_commands": [
            "formal-sensitivity --samples 1024 --workers {} --resume".format(workers),
            "design-matrix --workers {} --resume".format(workers),
            "counterfactual-bottleneck --workers {} --resume".format(workers),
            "design-experiments --method numerical-d-optimal",
        ],
        "evidence_label": UNCALIBRATED,
    }
    (base / "release_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    write_manifest(base / "release_manifest.json", finish_manifest(manifest, status="initialized"))
    return status

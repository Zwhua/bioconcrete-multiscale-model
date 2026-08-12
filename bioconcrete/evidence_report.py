"""Generate an evidence-status report without promoting missing results."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict


def evidence_report(project_root: Path, run_dir: Path, output_dir: Path) -> Dict[str, object]:
    public_root = project_root / "data" / "public"
    paths = {
        "public_calibration": run_dir / "frozen_run.json",
        "external_validation": project_root / "model_runs" / "external_validation" / "external_validation.json",
        "measurement_error": project_root / "model_runs" / "measurement_error" / "measurement_error.json",
        "formal_sensitivity": project_root / "model_runs" / "formal_sensitivity" / "formal_sensitivity.json",
        "design_matrix": project_root / "model_runs" / "design_matrix" / "design_summary.json",
    }
    status = {name: ("complete" if path.exists() else "not_available") for name, path in paths.items()}
    receipts = sorted(path.stem for path in (public_root / "receipts").glob("*.json")) if (public_root / "receipts").exists() else []
    report = {
        "project_experiment": {"rows": 0, "status": "not_available"},
        "public_data_receipts": receipts,
        "evidence_components": status,
        "model_prediction_label": "public-data-supported prediction requiring prospective validation",
        "claims": {
            "team_wet_lab_data_used": False,
            "virtual_dbtl_claimed": False,
            "external_parameters_refitted": False,
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "evidence_status.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    lines = ["# Evidence report", "", "The team currently has no project wet-lab time series.", ""]
    labels = {
        "public_calibration": "Public calibration data",
        "external_validation": "Independent external validation",
        "measurement_error": "Crack-width measurement error",
        "formal_sensitivity": "Formal SALib sensitivity",
        "design_matrix": "Preregistered design matrix",
    }
    for name, value in status.items():
        lines.append("- {}: `{}`".format(labels[name], value))
    lines.extend(["", "Missing components are not replaced with synthetic scientific evidence.",
                  "Synthetic tests validate software behavior only."])
    (output_dir / "EVIDENCE_REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report

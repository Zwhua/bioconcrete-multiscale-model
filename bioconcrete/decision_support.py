"""Generate evidence-bounded model-informed experimental decision tables."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict

import pandas as pd

REQUIRED_FIELDS = (
    "recommended_condition", "control_condition", "predicted_benefit",
    "uncertainty", "major_risk", "applicability", "required_wet_lab",
    "success_threshold", "failure_threshold", "evidence_level",
    "rationale", "config_hash", "code_hash",
)


def generate_decision_support(
    design_matrix_path: Path, output_dir: Path, config_hash: str, code_hash: str,
    counterfactual_path: Path = None,
) -> Dict[str, object]:
    """Convert model scenarios into recommendations without validation claims."""

    frame = pd.read_csv(design_matrix_path)
    if frame.empty:
        raise ValueError("Design matrix is empty")
    counterfactual = pd.read_csv(counterfactual_path) if counterfactual_path and counterfactual_path.exists() else pd.DataFrame()
    dominant = (
        str(counterfactual.sort_values("aggregate_control", ascending=False).iloc[0]["factor"])
        if len(counterfactual) else "pending_counterfactual_analysis"
    )
    records = []
    for _, row in frame.iterrows():
        records.append({
            "scenario_id": hashlib.sha256(row.to_json().encode()).hexdigest()[:12],
            "recommended_condition": "width={};dose={};wet={}h;activity={};delay={}h;leak={}".format(
                row.get("crack_width_mm"), row.get("agent_dosage"), row.get("wet_hours_per_day"),
                row.get("activity_multiplier"), row.get("response_delay_h"), row.get("basal_leak_fraction")),
            "control_condition": "default model configuration",
            "predicted_benefit": row.get("closure_28d"),
            "uncertainty": "prior predictive interval required",
            "major_risk": dominant,
            "applicability": "0.1-0.5 mm cracks; configured prior domain",
            "required_wet_lab": "closure plus CaCO3 mass and substrate at matched time points",
            "success_threshold": 0.5, "failure_threshold": 0.1,
            "evidence_level": "uncalibrated prospective model prediction",
            "rationale": "Pareto rank with {} counterfactual status".format(dominant),
            "config_hash": config_hash, "code_hash": code_hash,
            "decision": row.get("decision"), "dominant_bottleneck": dominant,
        })
    decisions = pd.DataFrame(records)
    output_dir.mkdir(parents=True, exist_ok=True)
    decisions.loc[decisions["decision"] == "recommended"].to_csv(output_dir / "recommended_conditions.csv", index=False)
    decisions.loc[decisions["decision"] == "not_recommended"].to_csv(output_dir / "rejected_conditions.csv", index=False)
    decisions.to_csv(output_dir / "tradeoff_table.csv", index=False)
    priorities = pd.DataFrame([
        {"measurement": "paired CaCO3 mass and closure", "priority": 1, "reason": "identify wall deposition and avoid geometry confounding"},
        {"measurement": "substrate concentration time series", "priority": 2, "reason": "identify release versus effective activity"},
        {"measurement": "dissolved oxygen time series", "priority": 3, "reason": "test oxygen/continuous-wetting assumption"},
        {"measurement": "calcium and DIC", "priority": 4, "reason": "separate calcium/carbon limitation"},
        {"measurement": "abiotic C-S-H closure control", "priority": 5, "reason": "separate initial filling from biomineralization"},
    ])
    priorities.to_csv(output_dir / "measurement_priorities.csv", index=False)
    (output_dir / "decision_summary.md").write_text(
        "# Model-informed experimental plan\n\nProspective DBTL loop. Awaiting wet-lab execution. Not experimentally validated.\n",
        encoding="utf-8",
    )
    provenance = {"config_hash": config_hash, "code_hash": code_hash,
                  "evidence_class": "uncalibrated prospective model prediction", "required_fields": list(REQUIRED_FIELDS)}
    (output_dir / "decision_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return {"scenarios": len(decisions), "recommended": int((decisions["decision"] == "recommended").sum())}

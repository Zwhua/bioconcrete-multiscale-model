"""Anonymous design-category to measurable model-parameter mapping."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import pandas as pd

from .evidence_state import UNCALIBRATED


MAPPINGS: List[Dict[str, object]] = [
    {"design_category": "surface_localization", "model_parameter": "kinetics.active_unit_concentration",
     "unit": "relative active unit m-3", "expected_direction": "increase",
     "required_measurement": "surface-associated signal and retained activity",
     "source_class": "design hypothesis", "identifiability": "weakly identifiable",
     "distinguishable_now": False},
    {"design_category": "linker_accessibility", "model_parameter": "kinetics.activity_multiplier",
     "unit": "dimensionless", "expected_direction": "context dependent",
     "required_measurement": "matched construct relative activity under concrete-relevant conditions",
     "source_class": "design hypothesis", "identifiability": "not identifiable without matched controls",
     "distinguishable_now": False},
    {"design_category": "candidate_enzyme_activity", "model_parameter": "kinetics.effective_kcat_s",
     "unit": "s-1", "expected_direction": "increase",
     "required_measurement": "aggregate kinetics and stability under matched pH and temperature",
     "source_class": "database prior", "identifiability": "scenario only",
     "distinguishable_now": False},
    {"design_category": "repair_unit_loading", "model_parameter": "kinetics.spore_density_rel",
     "unit": "relative inventory m-3", "expected_direction": "increase",
     "required_measurement": "viable or culturable unit inventory after encapsulation",
     "source_class": "project measurement required", "identifiability": "weakly identifiable",
     "distinguishable_now": False},
    {"design_category": "encapsulation_material", "model_parameter": "kinetics.capsule_release_s",
     "unit": "s-1", "expected_direction": "design dependent",
     "required_measurement": "release and retained-activity time courses",
     "source_class": "public calibration pending", "identifiability": "estimable with substrate time series",
     "distinguishable_now": True},
    {"design_category": "encapsulation_material", "model_parameter": "kinetics.basal_leak_fraction",
     "unit": "dimensionless", "expected_direction": "decrease",
     "required_measurement": "dry-state and pre-crack material loss",
     "source_class": "public calibration pending", "identifiability": "estimable with early inventory",
     "distinguishable_now": True},
    {"design_category": "csh_payload", "model_parameter": "kinetics.capsule_csh_volume_fraction",
     "unit": "m3 m-3 repair inventory", "expected_direction": "increase",
     "required_measurement": "abiotic C-S-H control and released solid volume",
     "source_class": "project measurement required", "identifiability": "confounded with wall deposition",
     "distinguishable_now": True},
]


def generate_biological_design(output_dir: Path) -> Dict[str, object]:
    """Write auditable anonymous mappings and prospective candidate decisions."""

    output_dir.mkdir(parents=True, exist_ok=True)
    mapping = pd.DataFrame(MAPPINGS)
    mapping["evidence_label"] = UNCALIBRATED
    mapping.to_csv(output_dir / "part_parameter_mapping.csv", index=False)
    candidates = pd.DataFrame([
        {"candidate_id": "anonymous_activity_baseline", "design_category": "candidate_enzyme_activity",
         "parameter_multiplier": 1.0, "candidate_status": "prior scoring only"},
        {"candidate_id": "anonymous_activity_high", "design_category": "candidate_enzyme_activity",
         "parameter_multiplier": 2.0, "candidate_status": "prior scoring only"},
        {"candidate_id": "low_leak_encapsulation", "design_category": "encapsulation_material",
         "parameter_multiplier": 0.5, "candidate_status": "prospective material class"},
        {"candidate_id": "no_csh_control", "design_category": "csh_payload",
         "parameter_multiplier": 0.0, "candidate_status": "required control"},
    ])
    candidates["evidence_label"] = UNCALIBRATED
    candidates.to_csv(output_dir / "candidate_constructs.csv", index=False)
    predictions = candidates.copy()
    predictions["prediction_type"] = "scenario parameter multiplier; not measured performance"
    predictions["requires_measurement"] = True
    predictions.to_csv(output_dir / "construct_predictions.csv", index=False)
    (output_dir / "biological_design_decisions.md").write_text(
        "# Anonymous biological/material design decisions\n\n"
        "Evidence level: uncalibrated prospective model prediction.\n\n"
        "The current model does not distinguish candidate sequences or construction details. "
        "It prioritizes encapsulation release/leak and C-S-H controls over broad activity screening. "
        "Any activity multiplier is a scenario prior until measured under matched conditions.\n",
        encoding="utf-8",
    )
    summary = {"mapping_rows": len(mapping), "candidate_rows": len(candidates),
               "contains_sequences": False, "evidence_label": UNCALIBRATED}
    (output_dir / "biological_design_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

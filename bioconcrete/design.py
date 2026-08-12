"""Preregistered prospective design matrix and transparent Pareto ranking."""

from __future__ import annotations

import copy
import itertools
import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import yaml

from .config import ModelConfig
from .model import simulate_0d


def _dominates(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.all(a <= b) and np.any(a < b))


def design_matrix(preregister_path: Path, output_dir: Path,
                  config: Optional[ModelConfig] = None, limit: Optional[int] = None) -> Dict[str, object]:
    prereg = yaml.safe_load(preregister_path.read_text(encoding="utf-8"))
    if prereg.get("status") != "fixed_before_external_validation":
        raise ValueError("Scenario file is not marked fixed before external validation")
    factors = prereg["factors"]
    names = list(factors)
    combinations = list(itertools.product(*(factors[name] for name in names)))
    if limit is not None:
        combinations = combinations[:max(int(limit), 0)]
    rows = []
    for values in combinations:
        setting = dict(zip(names, values))
        trial = copy.deepcopy(config or ModelConfig())
        trial.transport.crack_width_mm = setting["crack_width_mm"]
        trial.environment.wet_hours_per_day = setting["wet_hours_per_day"]
        trial.kinetics.activity_multiplier = setting["activity_multiplier"]
        trial.kinetics.response_delay_h = setting["response_delay_h"]
        trial.kinetics.basal_leak_fraction = setting["basal_leak_fraction"]
        trial.kinetics.agent_dosage_multiplier = setting["agent_dosage"]
        trial.simulation.days = 28.0
        trial.simulation.output_interval_days = 1.0
        result = simulate_0d(trial)
        frame = result.frame
        closure = frame["crack_closure_ratio"].to_numpy(float)
        reached = frame.loc[closure >= 0.5, "time_d"]
        premature = float(frame.loc[frame["time_d"] <= 1.0, "capsule_calcium_lactate_mol_m3"].iloc[-1])
        initial_inventory = float(frame["capsule_calcium_lactate_mol_m3"].iloc[0])
        rows.append({**setting,
            "closure_28d": result.summary["mean_crack_closure_ratio"],
            "time_to_50pct_d": float(reached.iloc[0]) if len(reached) else np.nan,
            "permeability_ratio": result.summary["mean_permeability_ratio"],
            "closure_per_agent": result.summary["mean_crack_closure_ratio"] / setting["agent_dosage"],
            "premature_consumption": 1.0 - premature / max(initial_inventory, 1e-30),
            "target_probability": np.nan,
            "evidence_label": prereg["evidence_label"],
        })
    frame = pd.DataFrame(rows)
    if len(frame):
        objectives = np.column_stack([
            -frame["closure_28d"], frame["permeability_ratio"],
            -frame["closure_per_agent"], frame["premature_consumption"],
        ])
        dominated = np.array([any(_dominates(objectives[j], objectives[i]) for j in range(len(frame)) if i != j)
                              for i in range(len(frame))])
        frame["decision"] = np.where(~dominated, "recommended",
                                     np.where(frame["closure_28d"] >= frame["closure_28d"].median(),
                                              "robust_alternative", "not_recommended"))
    output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_dir / "design_matrix.csv", index=False)
    summary = {"scenario_count": len(frame), "preregister_sha256": __import__("hashlib").sha256(
        preregister_path.read_bytes()).hexdigest(), "evidence_label": prereg["evidence_label"],
        "monte_carlo_target_probability_available": False}
    (output_dir / "design_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

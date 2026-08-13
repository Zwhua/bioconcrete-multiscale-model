"""Preregistered, resumable and parallel prospective design matrix."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import itertools
import json
from pathlib import Path
import time
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
import yaml

from .config import ModelConfig
from .model import simulate_0d
from .evidence_state import validate_evidence_label


def _dominates(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.all(a <= b) and np.any(a < b))


def _scenario_id(setting: Dict[str, object]) -> str:
    payload = json.dumps(setting, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _config_from_dict(raw: Dict[str, object]) -> ModelConfig:
    config = ModelConfig()
    for section, values in raw.items():
        for name, value in values.items():
            setattr(getattr(config, section), name, value)
    config.validate()
    return config


def _evaluate_scenario(payload: Tuple[Dict[str, object], Dict[str, object]]) -> Dict[str, object]:
    setting, raw_config = payload
    trial = _config_from_dict(raw_config)
    trial.transport.crack_width_mm = setting["crack_width_mm"]
    trial.environment.wet_hours_per_day = setting["wet_hours_per_day"]
    trial.kinetics.activity_multiplier = setting["activity_multiplier"]
    trial.kinetics.response_delay_h = setting["response_delay_h"]
    trial.kinetics.basal_leak_fraction = setting["basal_leak_fraction"]
    trial.kinetics.agent_dosage_multiplier = setting["agent_dosage"]
    trial.simulation.days = 28.0
    trial.simulation.output_interval_days = 1.0
    started = time.perf_counter()
    result = simulate_0d(trial)
    frame = result.frame
    closure = frame["crack_closure_ratio"].to_numpy(float)
    reached = frame.loc[closure >= 0.5, "time_d"]
    first_day = frame.loc[frame["time_d"] <= 1.0]
    initial = float(frame["capsule_calcium_lactate_mol_m3"].iloc[0])
    remaining = float(first_day["capsule_calcium_lactate_mol_m3"].iloc[-1])
    row = {
        **setting, "scenario_id": _scenario_id(setting),
        "closure_28d": result.summary["mean_crack_closure_ratio"],
        "time_to_50pct_d": float(reached.iloc[0]) if len(reached) else np.nan,
        "permeability_ratio": result.summary["mean_permeability_ratio"],
        "closure_per_agent": result.summary["mean_crack_closure_ratio"] / setting["agent_dosage"],
        "premature_consumption": 1.0 - remaining / max(initial, 1e-30),
        "target_probability": np.nan, "runtime_s": time.perf_counter() - started,
    }
    row["dominant_bottleneck"] = "pending_counterfactual_analysis"
    row["bottleneck_score"] = np.nan
    return row


def design_matrix(
    preregister_path: Path, output_dir: Path,
    config: Optional[ModelConfig] = None, limit: Optional[int] = None,
    workers: int = 1, resume: bool = False,
) -> Dict[str, object]:
    """Evaluate a fixed scenario matrix with deterministic IDs and resume."""

    prereg = yaml.safe_load(preregister_path.read_text(encoding="utf-8"))
    if prereg.get("status") != "fixed_before_external_validation":
        raise ValueError("Scenario file is not marked fixed before external validation")
    evidence_label = validate_evidence_label(prereg["evidence_label"])
    factors = prereg["factors"]
    names = list(factors)
    combinations = list(itertools.product(*(factors[name] for name in names)))
    if limit is not None:
        combinations = combinations[:max(int(limit), 0)]
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "design_matrix.csv"
    existing = pd.read_csv(output_path) if resume and output_path.exists() else pd.DataFrame()
    completed = set(existing.get("scenario_id", pd.Series(dtype=str)).astype(str))
    rows = existing.to_dict("records") if len(existing) else []
    raw_config = (config or ModelConfig()).to_dict()
    payloads = []
    for values in combinations:
        setting = dict(zip(names, values))
        if _scenario_id(setting) not in completed:
            payloads.append((setting, raw_config))
    failures = []

    def flush() -> None:
        table = pd.DataFrame(rows).sort_values("scenario_id")
        temporary = output_path.with_suffix(".tmp")
        table.to_csv(temporary, index=False)
        temporary.replace(output_path)
    pending = {"count": 0}

    def record(row: Dict[str, object]) -> None:
        rows.append(row)
        pending["count"] += 1
        if pending["count"] >= 16:
            flush()
            pending["count"] = 0

    if max(int(workers), 1) == 1:
        for payload in payloads:
            try:
                record(_evaluate_scenario(payload))
            except Exception as error:
                failures.append({"scenario_id": _scenario_id(payload[0]), "error": str(error)})
    else:
        with ProcessPoolExecutor(max_workers=int(workers)) as executor:
            futures = {executor.submit(_evaluate_scenario, payload): payload for payload in payloads}
            for future in as_completed(futures):
                try:
                    record(future.result())
                except Exception as error:
                    failures.append({"scenario_id": _scenario_id(futures[future][0]), "error": str(error)})

    if pending["count"]:
        flush()
    frame = pd.DataFrame(rows)
    if len(frame):
        frame = frame.sort_values("scenario_id").drop_duplicates("scenario_id", keep="last").reset_index(drop=True)
        frame["evidence_label"] = evidence_label
        objectives = np.column_stack([
            -frame["closure_28d"], frame["permeability_ratio"],
            -frame["closure_per_agent"], frame["premature_consumption"],
        ])
        dominated = np.array([
            any(_dominates(objectives[j], objectives[i]) for j in range(len(frame)) if i != j)
            for i in range(len(frame))
        ])
        frame["decision"] = np.where(
            ~dominated, "recommended",
            np.where(frame["closure_28d"] >= frame["closure_28d"].median(),
                     "robust_alternative", "not_recommended"),
        )
    frame.to_csv(output_path, index=False)
    pd.DataFrame(failures, columns=["scenario_id", "error"]).to_csv(
        output_dir / "failed_scenarios.csv", index=False
    )
    summary = {
        "scenario_count": len(frame),
        "preregister_sha256": hashlib.sha256(preregister_path.read_bytes()).hexdigest(),
        "evidence_label": evidence_label,
        "monte_carlo_target_probability_available": False,
        "workers": int(workers), "resume": bool(resume),
        "failed_scenarios": len(failures),
    }
    (output_dir / "design_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

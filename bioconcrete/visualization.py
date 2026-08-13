"""Reproducible static scientific figures from completed V5 artifacts only."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .evidence_state import UNCALIBRATED


COLORS = {"calcite": "#0072B2", "csh": "#E69F00", "other": "#009E73", "missing": "#999999"}


def _save(fig: plt.Figure, base: Path) -> None:
    fig.tight_layout()
    fig.savefig(base.with_suffix(".png"), dpi=220)
    fig.savefig(base.with_suffix(".svg"))
    plt.close(fig)


def _provenance(run_dir: Path) -> str:
    manifest = run_dir / "release_manifest.json"
    if not manifest.exists():
        return "v0.5.0 | manifest missing | {}".format(UNCALIBRATED)
    data = json.loads(manifest.read_text(encoding="utf-8"))
    return "v{} | {} | {} | {}".format(
        data.get("model_version", "0.5.0"), str(data.get("git_commit", "unknown"))[:8],
        str(data.get("config_sha256", "unknown"))[:10], data.get("evidence_label", UNCALIBRATED),
    )


def _decision_flow(output: Path, footer: str) -> None:
    labels = ["Design categories", "Measurable parameters", "Reaction/transport",
              "Mineral deposition", "Crack closure", "Counterfactual control", "Next experiment"]
    fig, axis = plt.subplots(figsize=(12, 2.6)); axis.set_axis_off()
    for index, label in enumerate(labels):
        x = index / (len(labels) - 1)
        axis.text(x, .55, label, ha="center", va="center", fontsize=9,
                  bbox={"boxstyle": "round,pad=.35", "fc": "#E8F1ED", "ec": "#24745E"})
        if index < len(labels) - 1:
            axis.annotate("", (x + 0.07, .55), (x + 0.02, .55), arrowprops={"arrowstyle": "->"})
    axis.set_title("Model-to-decision evidence chain")
    axis.text(0, .05, footer, fontsize=7, transform=axis.transAxes)
    _save(fig, output)


def _waterfall(summary_path: Path, output: Path, footer: str) -> None:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    csh = 100 * float(summary["mean_csh_closure_contribution"])
    calcite = 100 * float(summary["mean_calcite_closure_contribution"])
    total = 100 * float(summary["mean_crack_closure_ratio"])
    data = pd.DataFrame({"component": ["C-S-H", "CaCO3", "Unclosed"],
                         "percentage_points": [csh, calcite, max(100-total, 0)]})
    data.to_csv(output.with_name(output.name + "_data.csv"), index=False)
    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.bar(data["component"], data["percentage_points"], color=[COLORS["csh"], COLORS["calcite"], "#BBBBBB"])
    axis.set(ylabel="Share of initial crack width (percentage points)",
             title="Default closure is dominated by the C-S-H prior")
    axis.text(0, -.20, footer, fontsize=7, transform=axis.transAxes)
    _save(fig, output)


def _tornado(path: Path, output: Path, footer: str) -> None:
    frame = pd.read_csv(path).sort_values("aggregate_control")
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.barh(frame["factor"], frame["aggregate_control"], color="#0072B2")
    axis.set(xlabel="Median absolute normalized control coefficient", title="Counterfactual model controls")
    axis.text(0, -.18, footer, fontsize=7, transform=axis.transAxes)
    _save(fig, output)


def _uncertainty(path: Path, output: Path, footer: str) -> None:
    frame = pd.read_csv(path)
    y = np.arange(len(frame)); median = frame["median"].to_numpy(float)
    low, high = frame["lower_95"].to_numpy(float), frame["upper_95"].to_numpy(float)
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.errorbar(median, y, xerr=[median-low, high-median], fmt="o", color="#0072B2")
    axis.set(yticks=y, yticklabels=frame["metric"], xlabel="Model output (metric-specific units)",
             title="Prior predictive intervals: uncalibrated, not experimental data")
    axis.text(0, -.18, footer, fontsize=7, transform=axis.transAxes)
    _save(fig, output)


def _response_surface(path: Path, output: Path, footer: str) -> None:
    frame = pd.read_csv(path)
    selected = frame.groupby(["crack_width_mm", "agent_dosage"], as_index=False)["closure_28d"].median()
    pivot = selected.pivot(index="agent_dosage", columns="crack_width_mm", values="closure_28d")
    fig, axis = plt.subplots(figsize=(7, 5)); image = axis.imshow(100*pivot, origin="lower", aspect="auto", cmap="viridis")
    axis.set(xticks=range(len(pivot.columns)), xticklabels=pivot.columns, yticks=range(len(pivot.index)),
             yticklabels=pivot.index, xlabel="Initial crack width (mm)", ylabel="Complete agent dose multiplier",
             title="Width-dose response surface")
    fig.colorbar(image, ax=axis, label="28-day closure (%)"); axis.text(0, -.18, footer, fontsize=7, transform=axis.transAxes)
    _save(fig, output)


def _pareto(path: Path, output: Path, footer: str) -> None:
    frame = pd.read_csv(path)
    fig, axis = plt.subplots(figsize=(7, 5))
    groups = {name: group for name, group in frame.groupby("decision")}
    for name, group in groups.items():
        axis.scatter(group["agent_dosage"], 100*group["closure_28d"], label=name, alpha=.7)
    axis.set(xlabel="Complete agent dose multiplier", ylabel="28-day closure (%)", title="Pareto design trade-offs")
    axis.legend(); axis.text(0, -.18, footer, fontsize=7, transform=axis.transAxes); _save(fig, output)


def _information_gain(path: Path, output: Path, footer: str) -> None:
    frame = pd.read_csv(path).sort_values("rank").head(10)
    fig, axis = plt.subplots(figsize=(8, 5)); axis.barh(frame["experiment_id"].astype(str), frame["information_gain"], color="#009E73")
    axis.invert_yaxis(); axis.set(xlabel="Incremental log-determinant gain", title="Recommended complementary experiments")
    axis.text(0, -.18, footer, fontsize=7, transform=axis.transAxes); _save(fig, output)


def _multiscale(run_dir: Path, output: Path, footer: str) -> bool:
    files = list((run_dir / "baseline").rglob("state.csv"))
    if not files:
        return False
    fig, axes = plt.subplots(1, len(files), figsize=(5*len(files), 4), squeeze=False)
    for axis, path in zip(axes[0], files):
        frame = pd.read_csv(path); level = "2D" if "y_mm" in frame else ("1D" if "x_mm" in frame else "0D")
        grouped = frame.groupby("time_d")["crack_closure_ratio"].mean()
        axis.plot(grouped.index, 100*grouped.values); axis.set(title=level, xlabel="Time (d)", ylabel="Mean closure (%)")
    axes[0][0].text(0, -.25, footer, fontsize=7, transform=axes[0][0].transAxes); _save(fig, output)
    return True


def render_figures(run_dir: Path, output_dir: Optional[Path] = None) -> Dict[str, object]:
    """Render all figures supported by existing artifacts and report missing inputs."""

    output = output_dir or run_dir / "figures"; output.mkdir(parents=True, exist_ok=True)
    footer = _provenance(run_dir); generated, missing = [], []
    _decision_flow(output / "figure01_model_to_decision", footer); generated.append("figure01_model_to_decision")
    tasks = [
        ("figure02_closure_waterfall", next(iter((run_dir / "baseline").rglob("summary.json")), None), _waterfall),
        ("figure03_counterfactual_tornado", run_dir / "counterfactual_bottleneck" / "dominant_bottlenecks.csv", _tornado),
        ("figure04_prior_uncertainty", run_dir / "uncertainty" / "prior_predictive_summary.csv", _uncertainty),
        ("figure05_width_dose_surface", run_dir / "design_matrix" / "design_matrix.csv", _response_surface),
        ("figure06_pareto_front", run_dir / "design_matrix" / "design_matrix.csv", _pareto),
        ("figure07_experiment_information", run_dir / "experiment_design" / "recommended_experiments.csv", _information_gain),
    ]
    for name, path, function in tasks:
        if path is not None and Path(path).exists():
            function(Path(path), output / name, footer); generated.append(name)
        else:
            missing.append(name)
    if _multiscale(run_dir, output / "figure08_multiscale", footer): generated.append("figure08_multiscale")
    else: missing.append("figure08_multiscale")
    result = {"generated": generated, "missing": missing, "evidence_label": UNCALIBRATED,
              "substitute_data_generated": False}
    (output / "figure_manifest.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result

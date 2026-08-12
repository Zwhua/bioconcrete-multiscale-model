"""Figures and compact run reports."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _plot_0d(frame: pd.DataFrame, output: Path) -> None:
    fig, left = plt.subplots(figsize=(8, 5))
    left.plot(frame["time_d"], 100.0 * frame["crack_closure_ratio"], color="#176b87", linewidth=2, label="Closure")
    left.set(xlabel="Time (d)", ylabel="Crack closure (%)", ylim=(0, 100))
    right = left.twinx()
    right.plot(frame["time_d"], 100.0 * frame["transmissivity_ratio"], color="#b44c43", linewidth=2, label="Transmissivity")
    right.set(ylabel="Relative transmissivity (%)", ylim=(0, 105))
    handles = left.get_lines() + right.get_lines()
    left.legend(handles, [line.get_label() for line in handles], loc="center right")
    left.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _plot_1d(frame: pd.DataFrame, output: Path) -> None:
    fig, axis = plt.subplots(figsize=(8, 5))
    for day in sorted(frame["time_d"].unique()):
        selected = frame[frame["time_d"] == day]
        axis.plot(selected["x_mm"], 100.0 * selected["crack_closure_ratio"], linewidth=2, label="{:g} d".format(day))
    axis.set(xlabel="Position along crack (mm)", ylabel="Local crack closure (%)", ylim=(0, 100))
    axis.legend(ncol=2)
    axis.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def _plot_2d(frame: pd.DataFrame, output: Path) -> None:
    final = frame[frame["time_d"] == frame["time_d"].max()]
    pivot = final.pivot(index="y_mm", columns="x_mm", values="crack_closure_ratio")
    fig, axis = plt.subplots(figsize=(9, 3.4))
    image = axis.imshow(
        100.0 * pivot.to_numpy(),
        origin="lower",
        aspect="auto",
        extent=[final["x_mm"].min(), final["x_mm"].max(), final["y_mm"].min(), final["y_mm"].max()],
        cmap="viridis",
        vmin=0,
        vmax=max(10.0, float(100.0 * final["crack_closure_ratio"].max())),
    )
    axis.set(xlabel="Position along crack (mm)", ylabel="Crack width (mm)")
    colorbar = fig.colorbar(image, ax=axis)
    colorbar.set_label("Crack closure (%)")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def generate_report(run_dir: Path) -> Path:
    frame = pd.read_csv(run_dir / "state.csv")
    summary = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    diagnostics = json.loads((run_dir / "diagnostics.json").read_text(encoding="utf-8"))
    if "y_mm" in frame.columns:
        level = "2d"
        _plot_2d(frame, run_dir / "healing_map.png")
        figure = "healing_map.png"
    elif "x_mm" in frame.columns:
        level = "1d"
        _plot_1d(frame, run_dir / "healing_profiles.png")
        figure = "healing_profiles.png"
    else:
        level = "0d"
        _plot_0d(frame, run_dir / "timecourse.png")
        figure = "timecourse.png"
    report = """# BioConcrete model report

- Model level: `{level}`
- Final mean crack closure: {healing:.3%}
- Final maximum crack closure: {maximum:.3%}
- Mean permeability ratio: {permeability:.4f}
- Mean crack transmissivity ratio: {transmissivity:.4f}
- Mean calcite: {calcite:.3f} kg/m3 crack volume
- Maximum ammonium: {ammonium:.3e} mol/m3
- Nonnegative states: `{nonnegative}`
- Ammonia-free invariant: `{ammonia_free}`

![Model result]({figure})

The 80% healing value is an evaluation target, not a fitted or hard-coded outcome.
Database-derived values are priors. Project-specific rates require calibration against repair experiments.
""".format(
        level=level,
        healing=summary["mean_crack_closure_ratio"],
        maximum=summary["max_crack_closure_ratio"],
        permeability=summary["mean_permeability_ratio"],
        transmissivity=summary["mean_transmissivity_ratio"],
        calcite=summary["calcite_kg_m3_mean"],
        ammonium=summary["ammonium_mol_m3_max"],
        nonnegative=diagnostics["nonnegative"],
        ammonia_free=diagnostics["ammonia_free"],
        figure=figure,
    )
    path = run_dir / "REPORT.md"
    path.write_text(report, encoding="utf-8")
    return path

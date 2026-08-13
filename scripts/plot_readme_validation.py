"""Build the compact Gate D evidence panel used by the repository README."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = {"blue": "#176B87", "green": "#24745E", "gold": "#E69F00", "red": "#B34045"}


def _percent(comparison: dict[str, dict[str, float]]) -> tuple[list[str], np.ndarray]:
    labels = [
        "Calcite",
        "Area closure",
        "Max local closure",
        "Open volume",
        "O$_2$ penetration",
    ]
    keys = [
        "calcite_mol",
        "area_weighted_closure",
        "maximum_local_closure",
        "open_volume_m3",
        "oxygen_penetration_depth_m",
    ]
    return labels, 100.0 * np.asarray([comparison[key]["relative_error"] for key in keys])


def render(report_path: Path, output_path: Path) -> None:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    grid_labels, grid_error = _percent(report["grid_comparison_medium_fine"])
    time_labels, time_error = _percent(report["time_comparison_half_quarter"])

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.4), constrained_layout=True)
    fig.patch.set_facecolor("#F7FAF9")
    fig.suptitle("Gate D numerical verification", fontsize=18, fontweight="bold", color="#173B45")

    inventories = report["conservation"]
    conservation = [
        inventories["carbon_mol"]["relative_error"],
        inventories["calcium_mol"]["relative_error"],
    ]
    axes[0].bar(["Carbon", "Calcium"], conservation, color=[COLORS["blue"], COLORS["green"]])
    axes[0].axhline(5e-3, color=COLORS["red"], linestyle="--", label="0.5% limit")
    axes[0].set_yscale("log")
    axes[0].set_ylim(1e-17, 2e-2)
    axes[0].set_ylabel("Relative closure error")
    axes[0].set_title("Closed-system conservation")
    axes[0].legend(frameon=False, loc="upper right")

    x = np.arange(len(grid_labels))
    axes[1].bar(x, grid_error, color=COLORS["blue"])
    axes[1].axhline(5.0, color=COLORS["red"], linestyle="--", label="5% limit")
    axes[1].set_xticks(x, grid_labels, rotation=28, ha="right")
    axes[1].set_ylabel("Medium–fine difference (%)")
    axes[1].set_title("3D grid convergence")
    axes[1].set_ylim(0, 5.5)
    axes[1].legend(frameon=False, loc="upper right")

    axes[2].bar(x, time_error, color=COLORS["green"])
    axes[2].axhline(5.0, color=COLORS["red"], linestyle="--", label="5% limit")
    axes[2].set_xticks(x, time_labels, rotation=28, ha="right")
    axes[2].set_ylabel("3 h–1.5 h difference (%)")
    axes[2].set_title("Time-step convergence")
    axes[2].set_ylim(0, 5.5)
    axes[2].legend(frameon=False, loc="upper right")

    for axis in axes:
        axis.set_facecolor("white")
        axis.grid(axis="y", alpha=0.18)
        axis.spines[["top", "right"]].set_visible(False)

    fig.text(
        0.5,
        -0.03,
        "v0.6.0-development · Uncalibrated model output · Not experimental data",
        ha="center",
        fontsize=10,
        color="#53666B",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    render(args.report, args.output)


if __name__ == "__main__":
    main()

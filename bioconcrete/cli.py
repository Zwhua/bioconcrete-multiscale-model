"""Command-line interface for data preparation and model execution."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
from typing import Optional, Sequence

from .analysis import calibrate, sensitivity
from .chemistry import GeochemLookup, build_geochem_grid
from .config import ModelConfig
from .data_pipeline import prepare_data
from .model import simulate_0d, simulate_1d, simulate_2d
from .report import generate_report
from .validation import run_validation


def _root() -> Path:
    return Path(__file__).resolve().parents[1]


def _config(path: Optional[str]) -> ModelConfig:
    return ModelConfig.load(Path(path)) if path else ModelConfig()


def _lookup(path: Optional[str]) -> Optional[GeochemLookup]:
    if path:
        return GeochemLookup.load(Path(path))
    default = _root() / "data" / "processed" / "geochem" / "carbonate_lookup.csv"
    return GeochemLookup.load(default) if default.exists() else None


def _run_dir(base: Path, level: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return base / "{}_{}".format(stamp, level)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bioconcrete", description="Multiscale self-healing concrete model")
    subparsers = parser.add_subparsers(dest="command", required=True)

    config_parser = subparsers.add_parser("default-config", help="write a documented default JSON configuration")
    config_parser.add_argument("--output", default="model_config.json")

    prepare = subparsers.add_parser("prepare-data", help="build aggregate kinetic priors without biological design data")
    prepare.add_argument("--output", default="data/processed/model_priors")
    prepare.add_argument("--config")

    geochem = subparsers.add_parser("build-geochem-grid", help="build carbonate and cement-phase lookup files")
    geochem.add_argument("--output", default="data/processed/geochem")
    geochem.add_argument("--config")

    simulate = subparsers.add_parser("simulate", help="run the 0D, 1D, or true 2D model")
    simulate.add_argument("--level", choices=("0d", "1d", "2d"), required=True)
    simulate.add_argument("--config")
    simulate.add_argument("--geochem-grid")
    simulate.add_argument("--output", default="model_runs")

    calibration = subparsers.add_parser("calibrate", help="fit population, precipitation, and release rates")
    calibration.add_argument("--data", required=True)
    calibration.add_argument("--config")
    calibration.add_argument("--output", default="model_runs/calibration")
    calibration.add_argument("--bootstrap", type=int, default=20)

    sens = subparsers.add_parser("sensitivity", help="calculate Morris and Sobol indices")
    sens.add_argument("--config")
    sens.add_argument("--output", default="model_runs/sensitivity")
    sens.add_argument("--samples", type=int, default=8)

    report = subparsers.add_parser("report", help="create figures and a Markdown report for a simulation run")
    report.add_argument("--run", required=True)

    validation = subparsers.add_parser("validate", help="run conservation and discretization checks")
    validation.add_argument("--config")
    validation.add_argument("--output", default="model_runs/validation")
    validation.add_argument("--full", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    root = _root()
    if args.command == "default-config":
        path = Path(args.output)
        ModelConfig().save(path)
        print(path.resolve())
        return
    if args.command == "prepare-data":
        summary = prepare_data(root, Path(args.output), _config(args.config))
        print(json.dumps(summary, indent=2))
        return
    if args.command == "build-geochem-grid":
        config = _config(args.config)
        paths = build_geochem_grid(Path(args.output), config.chemistry, root, config.environment.temperature_c)
        print("\n".join(str(path.resolve()) for path in paths))
        return
    if args.command == "simulate":
        config = _config(args.config)
        lookup = _lookup(args.geochem_grid)
        runner = {"0d": simulate_0d, "1d": simulate_1d, "2d": simulate_2d}[args.level]
        result = runner(config, lookup)
        run_dir = _run_dir(Path(args.output), args.level)
        result.save(run_dir)
        generate_report(run_dir)
        print(run_dir.resolve())
        print(json.dumps(result.summary, indent=2))
        return
    if args.command == "calibrate":
        summary = calibrate(Path(args.data), Path(args.output), _config(args.config), bootstrap_samples=args.bootstrap)
        print(json.dumps(summary, indent=2))
        return
    if args.command == "sensitivity":
        morris, sobol = sensitivity(Path(args.output), _config(args.config), args.samples)
        print(sobol.sort_values("ST", ascending=False).to_string(index=False))
        return
    if args.command == "report":
        print(generate_report(Path(args.run)).resolve())
        return
    if args.command == "validate":
        result = run_validation(Path(args.output), _config(args.config), args.full)
        print(json.dumps(result, indent=2))

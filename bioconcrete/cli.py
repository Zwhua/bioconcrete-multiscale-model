"""Command-line interface for data preparation and model execution."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
from typing import Optional, Sequence

from .analysis import calibrate, sensitivity
from .chemistry import GeochemLookup, build_geochem_grid, compare_geochem_backends
from .config import ModelConfig
from .data_pipeline import prepare_data
from .model import simulate_0d, simulate_1d, simulate_2d
from .report import generate_report
from .validation import run_validation
from .evidence import calibrate_public, fit_measurement_error, validate_external
from .formal_analysis import formal_sensitivity
from .public_data import fetch_public_data, prepare_public_data
from .data_pipeline import parameter_registry
from .design import design_matrix
from .evidence_report import evidence_report
from .identifiability import identifiability_analysis
from .uncertainty import prior_predictive


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

    fetch = subparsers.add_parser("fetch-public-data", help="download manifest-selected public evidence files")
    fetch.add_argument("--manifest", default="data/public/DATASETS.yml")
    fetch.add_argument("--dataset")

    prepare_public = subparsers.add_parser("prepare-public-data", help="normalize one downloaded public dataset")
    prepare_public.add_argument("--dataset", required=True)

    geochem = subparsers.add_parser("build-geochem-grid", help="build carbonate and cement-phase lookup files")
    geochem.add_argument("--output", default="data/processed/geochem")
    geochem.add_argument("--config")
    phreeqc = subparsers.add_parser("build-phreeqc-grid", help="build a grid and record whether PHREEQC is actually available")
    phreeqc.add_argument("--output", default="data/processed/geochem")
    phreeqc.add_argument("--config")
    compare = subparsers.add_parser("compare-geochem-backends", help="compare PHREEQC and analytical backends when available")
    compare.add_argument("--grid", default="data/processed/geochem/carbonate_lookup.csv")
    compare.add_argument("--metadata", default="data/processed/geochem/geochem_metadata.json")
    compare.add_argument("--output", default="model_runs/geochem_comparison")

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

    formal = subparsers.add_parser("formal-sensitivity", help="run SALib Morris and Sobol analyses")
    formal.add_argument("--config")
    formal.add_argument("--output", default="model_runs/formal_sensitivity")
    formal.add_argument("--samples", type=int, default=256)

    identifiability = subparsers.add_parser(
        "identifiability", help="run practical local/Fisher identifiability diagnostics"
    )
    identifiability.add_argument("--config")
    identifiability.add_argument("--output", default="model_runs/identifiability")

    uncertainty = subparsers.add_parser(
        "prior-predictive", help="propagate parameter priors and scenario variability"
    )
    uncertainty.add_argument("--config")
    uncertainty.add_argument("--output", default="model_runs/prior_predictive")
    uncertainty.add_argument("--samples", type=int, default=256)
    uncertainty.add_argument("--seed", type=int, default=2026)
    uncertainty.add_argument("--resume", action="store_true")

    public_cal = subparsers.add_parser("calibrate-public", help="fit public calibration data with specimen holdout")
    public_cal.add_argument("--train", required=True)
    public_cal.add_argument("--config")
    public_cal.add_argument("--output", default="model_runs/public_calibration")
    public_cal.add_argument("--bootstrap", type=int, default=20)
    public_cal.add_argument("--profile-points", type=int, default=0)

    external = subparsers.add_parser("validate-external", help="validate a frozen public calibration run")
    external.add_argument("--dataset", required=True)
    external.add_argument("--frozen-run", required=True)
    external.add_argument("--output", default="model_runs/external_validation")

    measurement = subparsers.add_parser("fit-measurement-error", help="fit crack-width measurement noise only")
    measurement.add_argument("--dataset", required=True)
    measurement.add_argument("--output", default="model_runs/measurement_error")

    audit = subparsers.add_parser("audit-units", help="write the parameter/unit/source audit table")
    audit.add_argument("--config")
    audit.add_argument("--output", default="model_runs/unit_audit.csv")

    design = subparsers.add_parser("design-matrix", help="evaluate a preregistered prospective scenario matrix")
    design.add_argument("--preregister", default="PREREGISTERED_SCENARIOS.yml")
    design.add_argument("--config")
    design.add_argument("--output", default="model_runs/design_matrix")
    design.add_argument("--limit", type=int)

    evidence = subparsers.add_parser("evidence-report", help="report completed and missing evidence components")
    evidence.add_argument("--run", required=True)
    evidence.add_argument("--output", default="model_runs/evidence_report")

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
    if args.command == "fetch-public-data":
        summary = fetch_public_data(Path(args.manifest), root / "data" / "public", args.dataset)
        print(json.dumps(summary, indent=2))
        return
    if args.command == "prepare-public-data":
        summary = prepare_public_data(args.dataset, root / "data" / "public")
        print(json.dumps(summary, indent=2))
        return
    if args.command in {"build-geochem-grid", "build-phreeqc-grid"}:
        config = _config(args.config)
        paths = build_geochem_grid(Path(args.output), config.chemistry, root, config.environment.temperature_c)
        print("\n".join(str(path.resolve()) for path in paths))
        return
    if args.command == "compare-geochem-backends":
        result = compare_geochem_backends(Path(args.grid), Path(args.metadata), Path(args.output))
        print(json.dumps(result, indent=2))
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
    if args.command == "formal-sensitivity":
        _, sobol = formal_sensitivity(Path(args.output), _config(args.config), args.samples)
        print(sobol.sort_values("ST", ascending=False).to_string(index=False))
        return
    if args.command == "identifiability":
        result = identifiability_analysis(Path(args.output), _config(args.config))
        print(json.dumps(result, indent=2))
        return
    if args.command == "prior-predictive":
        result = prior_predictive(
            Path(args.output), _config(args.config), args.samples, args.seed, resume=args.resume
        )
        print(json.dumps(result, indent=2))
        return
    if args.command == "calibrate-public":
        result = calibrate_public(
            Path(args.train), Path(args.output), _config(args.config), args.bootstrap, args.profile_points
        )
        print(json.dumps(result, indent=2))
        return
    if args.command == "validate-external":
        result = validate_external(Path(args.dataset), Path(args.frozen_run), Path(args.output))
        print(json.dumps(result, indent=2))
        return
    if args.command == "fit-measurement-error":
        result = fit_measurement_error(Path(args.dataset), Path(args.output))
        print(json.dumps(result, indent=2))
        return
    if args.command == "audit-units":
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        parameter_registry(_config(args.config)).to_csv(path, index=False)
        print(path.resolve())
        return
    if args.command == "design-matrix":
        result = design_matrix(Path(args.preregister), Path(args.output), _config(args.config), args.limit)
        print(json.dumps(result, indent=2))
        return
    if args.command == "evidence-report":
        result = evidence_report(root, Path(args.run), Path(args.output))
        print(json.dumps(result, indent=2))
        return
    if args.command == "report":
        print(generate_report(Path(args.run)).resolve())
        return
    if args.command == "validate":
        result = run_validation(Path(args.output), _config(args.config), args.full)
        print(json.dumps(result, indent=2))

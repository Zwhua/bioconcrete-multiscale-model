"""Command-line interface for data preparation and model execution."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

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
from .experiment_design import rank_experiments
from .model_comparison import compare_structures
from .dashboard import generate_dashboard
from .decision_support import generate_decision_support
from .manifest import create_manifest, finish_manifest, write_manifest
from .counterfactual import counterfactual_bottleneck
from .biological_design import generate_biological_design
from .release_analysis import release_analysis
from .visualization import render_figures


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
    formal.add_argument("--workers", type=int, default=1)
    formal.add_argument("--resume", action="store_true")

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

    structures = subparsers.add_parser("compare-models", help="compare mechanistic and baseline structures")
    structures.add_argument("--config")
    structures.add_argument("--output", default="model_runs/model_comparison")
    structures.add_argument("--observations")

    experiments = subparsers.add_parser("design-experiments", help="rank prospective experiments by D-optimality")
    experiments.add_argument("--config")
    experiments.add_argument("--output", default="model_runs/experiment_design")
    experiments.add_argument("--method", choices=("numerical-d-optimal",), default="numerical-d-optimal")
    experiments.add_argument("--smoke", action="store_true")

    counterfactual = subparsers.add_parser(
        "counterfactual-bottleneck", help="calculate model-response control coefficients"
    )
    counterfactual.add_argument("--config")
    counterfactual.add_argument("--output", default="model_runs/counterfactual_bottleneck")
    counterfactual.add_argument("--perturbations", default="-0.2,-0.1,0.1,0.2")
    counterfactual.add_argument("--workers", type=int, default=1)
    counterfactual.add_argument("--resume", action="store_true")

    biological = subparsers.add_parser(
        "biological-design", help="export anonymous design-to-parameter mappings"
    )
    biological.add_argument("--output", default="model_runs/biological_design")

    release = subparsers.add_parser("release-analysis", help="initialize or resume v0.5.0 formal analyses")
    release.add_argument("--version", default="0.5.0")
    release.add_argument("--config")
    release.add_argument("--workers", type=int, default=16)
    release.add_argument("--resume", action="store_true")
    release.add_argument("--initialize-only", action="store_true")

    dashboard = subparsers.add_parser("dashboard", help="generate a read-only static evidence dashboard")
    dashboard.add_argument("--output", default="model_runs/dashboard")
    dashboard.add_argument("--run")

    figures = subparsers.add_parser("render-figures", help="render V5 scientific figures from completed artifacts")
    figures.add_argument("--run", default="model_runs/v0.5.0")
    figures.add_argument("--output")

    decisions = subparsers.add_parser("decision-support", help="build model-informed decision tables")
    decisions.add_argument("--design-matrix", default="model_runs/design_matrix/design_matrix.csv")
    decisions.add_argument("--config-hash", required=True)
    decisions.add_argument("--code-hash", required=True)
    decisions.add_argument("--output", default="model_runs/decision_support")
    decisions.add_argument("--counterfactual")

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
    design.add_argument("--workers", type=int, default=1)
    design.add_argument("--resume", action="store_true")

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
        manifest = create_manifest(root, _config(args.config), ["formal-sensitivity"], 2026)
        _, sobol = formal_sensitivity(
            Path(args.output), _config(args.config), args.samples, 2026, args.workers, args.resume
        )
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(sobol.sort_values("ST", ascending=False).to_string(index=False))
        return
    if args.command == "identifiability":
        manifest = create_manifest(root, _config(args.config), ["identifiability"], 2026)
        result = identifiability_analysis(Path(args.output), _config(args.config))
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(json.dumps(result, indent=2))
        return
    if args.command == "prior-predictive":
        manifest = create_manifest(root, _config(args.config), ["prior-predictive"], args.seed)
        result = prior_predictive(
            Path(args.output), _config(args.config), args.samples, args.seed, resume=args.resume
        )
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(json.dumps(result, indent=2))
        return
    if args.command == "compare-models":
        manifest = create_manifest(root, _config(args.config), ["compare-models"], 2026)
        observations = pd.read_csv(args.observations) if args.observations else None
        result = compare_structures(Path(args.output), _config(args.config), observations)
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(json.dumps(result, indent=2))
        return
    if args.command == "design-experiments":
        manifest = create_manifest(root, _config(args.config), ["design-experiments"], 2026)
        result = rank_experiments(Path(args.output), _config(args.config), smoke=args.smoke)
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(json.dumps(result, indent=2))
        return
    if args.command == "counterfactual-bottleneck":
        changes = tuple(float(value) for value in args.perturbations.split(","))
        manifest = create_manifest(root, _config(args.config), ["counterfactual-bottleneck"], 2026)
        result = counterfactual_bottleneck(
            Path(args.output), _config(args.config), changes, args.workers, args.resume
        )
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
        print(json.dumps(result, indent=2))
        return
    if args.command == "biological-design":
        result = generate_biological_design(Path(args.output))
        print(json.dumps(result, indent=2))
        return
    if args.command == "release-analysis":
        result = release_analysis(
            root, args.version, _config(args.config), args.workers, args.resume, args.initialize_only
        )
        print(json.dumps(result, indent=2))
        return
    if args.command == "dashboard":
        result = generate_dashboard(root, Path(args.output), Path(args.run) if args.run else None)
        print(json.dumps(result, indent=2))
        return
    if args.command == "render-figures":
        result = render_figures(Path(args.run), Path(args.output) if args.output else None)
        print(json.dumps(result, indent=2))
        return
    if args.command == "decision-support":
        result = generate_decision_support(
            Path(args.design_matrix), Path(args.output), args.config_hash, args.code_hash,
            Path(args.counterfactual) if args.counterfactual else None,
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
        manifest = create_manifest(
            root, _config(args.config), ["design-matrix"], 2026,
            {"preregister": Path(args.preregister)},
        )
        result = design_matrix(
            Path(args.preregister), Path(args.output), _config(args.config),
            args.limit, args.workers, args.resume,
        )
        write_manifest(Path(args.output) / "run_manifest.json", finish_manifest(manifest))
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

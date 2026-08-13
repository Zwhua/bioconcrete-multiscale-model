"""Interactive single-file evidence dashboard backed only by completed artifacts."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .evidence_state import UNCALIBRATED


def _table(path: Path, limit: int = 12) -> str:
    if not path.exists():
        return '<p class="missing">Missing evidence. No demonstration data were substituted.</p>'
    return pd.read_csv(path).head(limit).to_html(index=False, border=0, classes="data-table")


def _interactive_design(path: Path) -> str:
    """Return a self-contained Plotly view only when formal design data exist."""

    if not path.exists():
        return '<p class="missing">Interactive scenario view unavailable: design matrix missing.</p>'
    try:
        import plotly.express as px
        from plotly.offline import plot
    except ImportError:
        return '<p class="missing">Interactive view unavailable: install the pinned Plotly dependency.</p>'
    frame = pd.read_csv(path)
    required = {"agent_dosage", "closure_28d", "crack_width_mm", "wet_hours_per_day", "activity_multiplier"}
    if not required.issubset(frame):
        return '<p class="missing">Interactive view unavailable: required design fields are missing.</p>'
    figure = px.scatter(
        frame, x="agent_dosage", y="closure_28d", color="crack_width_mm",
        symbol="wet_hours_per_day", hover_data=["activity_multiplier", "scenario_id"],
        labels={"agent_dosage": "Complete agent dose multiplier", "closure_28d": "28-day closure ratio"},
        title="Scenario A/B and Pareto design explorer",
        color_continuous_scale="Viridis",
    )
    figure.update_layout(font={"family": "Arial"}, title_font={"size": 18})
    return plot(figure, include_plotlyjs=True, output_type="div", config={"responsive": True})


def generate_dashboard(project_root: Path, output_dir: Path, run_dir: Optional[Path] = None) -> Dict[str, object]:
    """Build a mobile-readable HTML dashboard and explicit static fallback."""

    run = run_dir or project_root / "model_runs" / "v0.5.0"
    assets = {
        "change": run / "biological_design" / "construct_predictions.csv",
        "controls": run / "counterfactual_bottleneck" / "dominant_bottlenecks.csv",
        "uncertainty": run / "uncertainty" / "prior_predictive_summary.csv",
        "design": run / "design_matrix" / "design_matrix.csv",
        "experiments": run / "experiment_design" / "recommended_experiments.csv",
        "missing": run / "release_status.json",
    }
    interactive = _interactive_design(assets["design"])
    questions = [
        ("What did the model change?", "change"), ("What controls healing?", "controls"),
        ("How uncertain are predictions?", "uncertainty"), ("Which design is recommended?", "design"),
        ("Which experiment should be run next?", "experiments"), ("What evidence is still missing?", "missing"),
    ]
    sections = []
    for title, key in questions:
        if key == "missing" and assets[key].exists():
            content = "<pre>{}</pre>".format(html.escape(assets[key].read_text(encoding="utf-8")))
        else:
            content = _table(assets[key])
        sections.append("<section id='{}'><h2>{}</h2>{}</section>".format(key, title, content))
    output_dir.mkdir(parents=True, exist_ok=True)
    document = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>BioConcrete v0.5.0 Evidence Dashboard</title>
<style>body{{font-family:Arial,sans-serif;margin:0;color:#17211d;background:#f7f8f6}}header{{background:#163d34;color:white;padding:24px 5vw}}
nav{{display:flex;gap:8px;overflow:auto;padding:10px 5vw;background:white;position:sticky;top:0}}nav a{{color:#174f42;white-space:nowrap}}
main{{max-width:1200px;margin:auto;padding:20px}}section{{padding:16px 0;border-bottom:1px solid #ccd6d1}}h1,h2{{letter-spacing:0}}
.evidence{{background:#e8f1ed;border-left:5px solid #24745e;padding:12px}}.missing{{color:#8a3f2c;font-weight:600}}
.data-table{{border-collapse:collapse;width:100%;background:white;display:block;overflow:auto}}th,td{{padding:7px;border:1px solid #d8ddd9;font-size:12px}}
pre{{white-space:pre-wrap;background:white;padding:12px}}@media(max-width:600px){{main{{padding:12px}}h1{{font-size:25px}}}}</style></head>
<body><header><h1>BioConcrete Model-to-Decision Dashboard</h1><p>Read-only v0.5.0 evidence system</p></header>
<nav>{}</nav><main><div class="evidence"><strong>Evidence:</strong> {}. Team wet-lab observations: 0. Public calibration: incomplete.</div>
<section><h2>Interactive design explorer</h2>{}</section>{}</main></body></html>""".format(
        "".join("<a href='#{}'>{}</a>".format(key, title) for title, key in questions),
        UNCALIBRATED, interactive, "".join(sections),
    )
    path = output_dir / "index.html"; path.write_text(document, encoding="utf-8")
    fallback = output_dir / "STATIC_FALLBACK.md"
    fallback.write_text("# BioConcrete evidence dashboard\n\nEvidence: {}.\n\n".format(UNCALIBRATED)
                        + "\n".join("- {}: {}".format(title, assets[key]) for title, key in questions), encoding="utf-8")
    manifest = {"dashboard": str(path), "static_fallback": str(fallback),
                "source_artifacts": {key: str(value) for key, value in assets.items()},
                "team_wet_lab_rows": 0, "generated_data": False,
                "evidence_label": UNCALIBRATED}
    (output_dir / "dashboard_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest

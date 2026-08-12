"""Generate a read-only static evidence dashboard from existing run artifacts."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Dict

import pandas as pd


def _table(path: Path, limit: int = 10) -> str:
    if not path.exists():
        return '<p class="missing">Not available. No substitute data were generated.</p>'
    frame = pd.read_csv(path).head(limit)
    return frame.to_html(index=False, border=0, classes="data-table")


def generate_dashboard(project_root: Path, output_dir: Path) -> Dict[str, object]:
    """Build a self-contained HTML dashboard without creating scientific data."""

    assets = {
        "design": project_root / "model_runs" / "design_matrix" / "design_matrix.csv",
        "uncertainty": project_root / "model_runs" / "prior_predictive" / "prior_predictive_summary.csv",
        "sensitivity": project_root / "model_runs" / "formal_sensitivity" / "sobol_indices.csv",
        "identifiability": project_root / "model_runs" / "identifiability" / "identifiability_table.csv",
        "experiments": project_root / "model_runs" / "experiment_design" / "recommended_experiments.csv",
        "applicability": project_root / "data" / "public" / "curated" / "applicability_domain.csv",
    }
    sections = []
    labels = {
        "design": "Design trade-offs and Pareto scenarios",
        "uncertainty": "Prior predictive uncertainty",
        "sensitivity": "Global sensitivity",
        "identifiability": "Practical identifiability",
        "experiments": "Recommended future experiments",
        "applicability": "Applicability domain",
    }
    for name, path in assets.items():
        sections.append("<section><h2>{}</h2>{}</section>".format(
            html.escape(labels[name]), _table(path)
        ))
    output_dir.mkdir(parents=True, exist_ok=True)
    document = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>BioConcrete Model Evidence Dashboard</title><style>
body{{font-family:Arial,sans-serif;margin:0;color:#17211d;background:#f7f8f6}}
header{{background:#163d34;color:white;padding:28px max(5vw,24px)}}
main{{max-width:1200px;margin:auto;padding:24px}} section{{margin:0 0 28px}}
h1,h2{{letter-spacing:0}} h2{{font-size:20px;border-bottom:2px solid #ca6b3f;padding-bottom:7px}}
.evidence{{background:#e8f1ed;border-left:5px solid #24745e;padding:14px;margin-top:18px}}
.missing{{color:#8a3f2c;font-weight:600}} .data-table{{border-collapse:collapse;width:100%;background:white}}
th,td{{padding:8px;border:1px solid #d8ddd9;text-align:left;font-size:13px}} th{{background:#e4ece8}}
</style></head><body><header><h1>BioConcrete Model Evidence Dashboard</h1>
<p>Read-only model outputs, priors, evidence gaps and future experiments</p></header><main>
<div class="evidence"><strong>Evidence boundary:</strong> Team wet-lab observations: 0.
Simulation tables are model output; parameter ranges are literature/scenario priors; missing public calibration remains missing.</div>
{}<section><h2>Mass conservation and model structure</h2><p>Run <code>python -m bioconcrete validate</code> and
<code>python -m bioconcrete compare-models</code> to regenerate current checks. Historical result images retain their producing version.</p></section>
</main></body></html>""".format("".join(sections))
    path = output_dir / "index.html"
    path.write_text(document, encoding="utf-8")
    manifest = {"dashboard": str(path), "source_artifacts": {k: str(v) for k, v in assets.items()},
                "team_wet_lab_rows": 0, "generated_data": False}
    (output_dir / "dashboard_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest

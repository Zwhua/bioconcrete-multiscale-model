"""Traceable acquisition and normalization of public concrete datasets."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Dict, Iterable, List, Mapping, Optional
import zipfile

import numpy as np
import pandas as pd
import requests
import yaml


OBSERVATION_COLUMNS = [
    "dataset_id", "specimen_id", "group_id", "split", "time_d",
    "initial_crack_width_mm", "current_crack_width_mm", "conditioning",
    "wet_hours_per_day", "agent_dosage", "calcite_mass_mg",
    "permeability_ratio", "stiffness_ratio", "measurement_sd",
    "source_file", "source_location",
]


def load_manifest(path: Path) -> Mapping[str, object]:
    content = yaml.safe_load(path.read_text(encoding="utf-8"))
    if content.get("schema_version") != 1 or not isinstance(content.get("datasets"), dict):
        raise ValueError("Unsupported public-data manifest")
    return content


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _record(record_id: int) -> Mapping[str, object]:
    response = requests.get(
        "https://zenodo.org/api/records/{}".format(record_id), timeout=60,
        headers={"User-Agent": "bioconcrete-evidence/0.2 (+https://github.com/Zwhua/bioconcrete-multiscale-model)"},
    )
    response.raise_for_status()
    return response.json()


def fetch_public_data(manifest_path: Path, root: Path, dataset_id: Optional[str] = None) -> Dict[str, object]:
    """Download only manifest-selected files and record checksums/licenses."""

    manifest = load_manifest(manifest_path)
    selected = manifest["datasets"]
    if dataset_id:
        if dataset_id not in selected:
            raise ValueError("Unknown dataset_id: {}".format(dataset_id))
        selected = {dataset_id: selected[dataset_id]}
    receipts = []
    for identifier, spec in selected.items():
        api_available = True
        try:
            metadata = _record(int(spec["record_id"]))
        except requests.HTTPError:
            api_available = False
            metadata = {"files": [], "metadata": {}}
        files = {item["key"]: item for item in metadata.get("files", [])}
        requested = list(spec.get("selected_files", []))
        if not requested:
            requested = [name for name in files if name.lower().endswith((".csv", ".xlsx"))]
        raw_dir = root / "raw" / identifier
        raw_dir.mkdir(parents=True, exist_ok=True)
        file_receipts = []
        for name in requested:
            item = files.get(name, {})
            destination = raw_dir / name
            links = item.get("links", {})
            url = links.get("self") or links.get("content") or (
                "https://zenodo.org/records/{}/files/{}?download=1".format(spec["record_id"], requests.utils.quote(name))
            )
            with requests.get(url, stream=True, timeout=120, headers={"User-Agent": "bioconcrete-evidence/0.2"}) as response:
                response.raise_for_status()
                with destination.open("wb") as handle:
                    for block in response.iter_content(1024 * 1024):
                        if block:
                            handle.write(block)
            file_receipts.append({
                "name": name,
                "url": url,
                "bytes": destination.stat().st_size,
                "sha256": _sha256(destination),
                "zenodo_checksum": item.get("checksum"),
            })
        receipt = {
            "dataset_id": identifier,
            "record_id": spec["record_id"],
            "doi": spec["doi"],
            "role": spec["role"],
            "license": metadata.get("metadata", {}).get("license", {}).get("id", "unspecified"),
            "zenodo_api_metadata_available": api_available,
            "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
            "files": file_receipts,
        }
        receipt_dir = root / "receipts"
        receipt_dir.mkdir(parents=True, exist_ok=True)
        (receipt_dir / "{}.json".format(identifier)).write_text(json.dumps(receipt, indent=2), encoding="utf-8")
        receipts.append(receipt)
    return {"datasets": receipts, "raw_data_committed": False}


def grouped_split(specimen_ids: Iterable[str], seed: int = 2026, holdout_fraction: float = 0.2) -> Dict[str, str]:
    """Stable split keyed only by specimen, preventing time-point leakage."""

    mapping = {}
    for specimen in sorted(set(str(item) for item in specimen_ids)):
        token = hashlib.sha256((str(seed) + ":" + specimen).encode("utf-8")).digest()
        fraction = int.from_bytes(token[:8], "big") / float(2**64)
        mapping[specimen] = "internal_test" if fraction < holdout_fraction else "train"
    return mapping


def _number(value: object) -> float:
    match = re.search(r"[-+]?\d*\.?\d+", str(value).replace(",", "."))
    return float(match.group()) if match else np.nan


def _blank(dataset_id: str, source_file: str, location: str) -> Dict[str, object]:
    row = {name: np.nan for name in OBSERVATION_COLUMNS}
    row.update(dataset_id=dataset_id, source_file=source_file, source_location=location)
    return row


def _spreadsheet_rows(path: Path, dataset_id: str) -> List[Dict[str, object]]:
    """Conservative generic extraction; ambiguous values stay missing."""

    rows: List[Dict[str, object]] = []
    workbook = pd.ExcelFile(path)
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet, header=None)
        for index, values in frame.iterrows():
            joined = " | ".join(str(value) for value in values if pd.notna(value))
            lower = joined.lower()
            if not joined or not any(token in lower for token in ("crack", "healing", "stiff", "mass")):
                continue
            row = _blank(dataset_id, path.name, "sheet={},row={}".format(sheet, index + 1))
            row["specimen_id"] = str(values.iloc[0]) if pd.notna(values.iloc[0]) else "{}:{}".format(sheet, index + 1)
            row["group_id"] = sheet
            numbers = [_number(value) for value in values]
            numbers = [value for value in numbers if np.isfinite(value)]
            if "crack" in lower and numbers:
                widths = [value for value in numbers if 0 < value <= 5]
                if widths:
                    row["current_crack_width_mm"] = widths[-1]
                    row["initial_crack_width_mm"] = widths[0]
            if "stiff" in lower and numbers:
                row["stiffness_ratio"] = numbers[-1] / 100.0 if numbers[-1] > 1 else numbers[-1]
            rows.append(row)
    return rows


def prepare_public_data(dataset_id: str, root: Path) -> Dict[str, object]:
    """Normalize available tabular files while retaining exact source locations."""

    raw_dir = root / "raw" / dataset_id
    if not raw_dir.exists():
        raise FileNotFoundError("Run fetch-public-data first: {}".format(dataset_id))
    work_dir = root / "work" / dataset_id
    local_files = sorted(path for path in raw_dir.rglob("*") if path.is_file())
    receipt_dir = root / "receipts"
    receipt_dir.mkdir(parents=True, exist_ok=True)
    local_receipt = {
        "dataset_id": dataset_id,
        "receipt_type": "local_manual_or_automated_files",
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": [{"name": str(path.relative_to(raw_dir)), "bytes": path.stat().st_size,
                   "sha256": _sha256(path)} for path in local_files],
        "publisher_checksum_verified": False,
    }
    (receipt_dir / "{}_local.json".format(dataset_id)).write_text(
        json.dumps(local_receipt, indent=2), encoding="utf-8"
    )
    if work_dir.exists():
        shutil.rmtree(str(work_dir))
    work_dir.mkdir(parents=True)
    for archive in raw_dir.glob("*.zip"):
        with zipfile.ZipFile(str(archive)) as handle:
            handle.extractall(str(work_dir))
    files = list(raw_dir.glob("*.xlsx")) + list(work_dir.rglob("*.xlsx"))
    rows: List[Dict[str, object]] = []
    for path in files:
        rows.extend(_spreadsheet_rows(path, dataset_id))
    frame = pd.DataFrame(rows, columns=OBSERVATION_COLUMNS)
    if not frame.empty:
        mapping = grouped_split(frame["specimen_id"].astype(str))
        frame["split"] = frame["specimen_id"].astype(str).map(mapping)
        if dataset_id == "marine_external":
            frame["split"] = "external_validation"
    output_dir = root / "derived" / dataset_id
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "observations.csv"
    frame.to_csv(output, index=False)
    summary = {
        "dataset_id": dataset_id,
        "rows": int(len(frame)),
        "files_scanned": [str(path.relative_to(root)) for path in files],
        "local_files_audited": int(len(local_files)),
        "local_receipt": str(receipt_dir / "{}_local.json".format(dataset_id)),
        "output": str(output),
        "warning": "Generic extraction retains only unambiguous fields; review source_location before calibration.",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

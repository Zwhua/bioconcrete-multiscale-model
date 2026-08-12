"""Prepare aggregate, non-sequence parameter priors from local public data."""

from __future__ import annotations

from dataclasses import asdict
import json
import mmap
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .config import ModelConfig, PARAMETER_PROVENANCE


TARGET_EC = ("1.1.1.27", "4.2.1.1")


def _quantile_rows(frame: pd.DataFrame, group_columns: Sequence[str], value_column: str, source: str) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    if frame.empty:
        return rows
    for keys, group in frame.groupby(list(group_columns), dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values = pd.to_numeric(group[value_column], errors="coerce").dropna()
        if values.empty:
            continue
        row = dict(zip(group_columns, keys))
        row.update(
            {
                "source": source,
                "count": int(values.size),
                "minimum": float(values.min()),
                "q25": float(values.quantile(0.25)),
                "median": float(values.median()),
                "q75": float(values.quantile(0.75)),
                "maximum": float(values.max()),
            }
        )
        rows.append(row)
    return rows


def summarize_sabio(raw_dir: Path) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for path in sorted(raw_dir.glob("*.json")):
        document = json.loads(path.read_text(encoding="utf-8"))
        records = document.get("data", [])
        if not records:
            continue
        ec = str(records[0].get("enzyme_description", {}).get("ec_number", ""))
        if ec not in TARGET_EC:
            continue
        parameters = []
        conditions = []
        for record in records:
            experiment = record.get("experimental_conditions", {})
            ph = (experiment.get("envvar_ph") or {}).get("start_value")
            temperature = (experiment.get("envvar_temperature") or {}).get("start_value")
            if ph is not None:
                conditions.append(("pH", float(ph)))
            if temperature is not None:
                conditions.append(("temperature_c", float(temperature)))
            for parameter in record.get("kineticlaw", {}).get("parameter", []) or []:
                value = parameter.get("n_start_value")
                if value is None:
                    continue
                parameter_type = (parameter.get("parameter_type") or {}).get("name") or parameter.get("name") or "unknown"
                unit = (parameter.get("unit") or {}).get("n_name") or (parameter.get("unit") or {}).get("name") or "unknown"
                parameters.append({"ec_number": ec, "parameter": parameter_type, "unit": unit, "value": value})
        rows.extend(_quantile_rows(pd.DataFrame(parameters), ["ec_number", "parameter", "unit"], "value", "SABIO-RK aggregate"))
        condition_frame = pd.DataFrame(conditions, columns=["parameter", "value"])
        if not condition_frame.empty:
            condition_frame["ec_number"] = ec
            condition_frame["unit"] = condition_frame["parameter"].map({"pH": "dimensionless", "temperature_c": "degree_C"})
            rows.extend(_quantile_rows(condition_frame, ["ec_number", "parameter", "unit"], "value", "SABIO-RK aggregate"))
    return pd.DataFrame(rows)


def _extract_top_level_json_value(path: Path, key: str) -> Optional[Any]:
    """Extract one large top-level enzyme object without loading the full BRENDA file."""

    token = ('"{}"'.format(key)).encode("utf-8")
    with path.open("rb") as handle:
        with mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as mapped:
            data_marker = mapped.find(b'"data"')
            position = mapped.find(token, max(data_marker, 0))
            while position >= 0:
                cursor = position + len(token)
                while cursor < len(mapped) and mapped[cursor] in b" \t\r\n":
                    cursor += 1
                if cursor < len(mapped) and mapped[cursor] == ord(":"):
                    break
                position = mapped.find(token, position + len(token))
            if position < 0:
                return None
            cursor += 1
            while mapped[cursor] in b" \t\r\n":
                cursor += 1
            opening = mapped[cursor]
            if opening not in (ord("{"), ord("[")):
                return None
            closing = ord("}") if opening == ord("{") else ord("]")
            depth = 0
            in_string = False
            escaped = False
            end = cursor
            while end < len(mapped):
                byte = mapped[end]
                if in_string:
                    if escaped:
                        escaped = False
                    elif byte == ord("\\"):
                        escaped = True
                    elif byte == ord('"'):
                        in_string = False
                elif byte == ord('"'):
                    in_string = True
                elif byte == opening:
                    depth += 1
                elif byte == closing:
                    depth -= 1
                    if depth == 0:
                        end += 1
                        break
                end += 1
            return json.loads(mapped[cursor:end].decode("utf-8"))


def _reported_numbers(value: Any) -> List[float]:
    if isinstance(value, (int, float)):
        return [float(value)]
    if not isinstance(value, str):
        return []
    prefix = value.split("{", 1)[0]
    return [float(item) for item in re.findall(r"(?<![A-Za-z])[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", prefix)]


def summarize_brenda(path: Path) -> pd.DataFrame:
    rows = []
    fields = {
        "km_value": "mM_as_reported",
        "turnover_number": "s^-1_as_reported",
        "kcat_km_value": "mM^-1_s^-1_as_reported",
        "ph_optimum": "dimensionless",
        "temperature_optimum": "degree_C",
    }
    for ec in TARGET_EC:
        enzyme = _extract_top_level_json_value(path, ec)
        if enzyme is None:
            continue
        values = []
        for field, unit in fields.items():
            for record in enzyme.get(field, []) or []:
                reported = record.get("value") if isinstance(record, dict) else record
                numbers = [number for number in _reported_numbers(reported) if number >= 0.0]
                if numbers:
                    # Ranges contribute their midpoint as one observation.
                    values.append((field, unit, float(np.mean(numbers[:2]))))
        if not values:
            rows.append({"ec_number": ec, "parameter": "records_present", "unit": "count", "source": "BRENDA aggregate", "count": 1, "minimum": np.nan, "q25": np.nan, "median": np.nan, "q75": np.nan, "maximum": np.nan})
            continue
        frame = pd.DataFrame(values, columns=["parameter", "unit", "value"])
        frame["ec_number"] = ec
        rows.extend(_quantile_rows(frame, ["ec_number", "parameter", "unit"], "value", "BRENDA aggregate"))
    return pd.DataFrame(rows)


def parameter_registry(config: ModelConfig) -> pd.DataFrame:
    values = config.to_dict()
    units = {
        "kinetics.qmax_lactate_mol_m3_s": "mol m^-3 s^-1",
        "kinetics.k_lactate_mol_m3": "mol m^-3",
        "kinetics.k_oxygen_mol_m3": "mol m^-3",
        "kinetics.decay_s": "s^-1",
        "kinetics.capsule_release_s": "s^-1",
        "chemistry.calcite_rate_mol_m3_s": "mol m^-3 s^-1",
        "chemistry.calcite_ksp": "mol^2 L^-2",
        "transport.diffusivity_oxygen_m2_s": "m^2 s^-1",
        "transport.crack_width_mm": "mm",
        "chemistry.portlandite_mol_m3": "mol m^-3",
    }
    rows = []
    for parameter, provenance in PARAMETER_PROVENANCE.items():
        section, name = parameter.split(".", 1)
        source_class, source_note, lower, upper = provenance
        rows.append(
            {
                "parameter": parameter,
                "default_value": values[section][name],
                "unit": units.get(parameter, "dimensionless"),
                "source_class": source_class,
                "source_note": source_note,
                "lower_bound": lower,
                "upper_bound": upper,
            }
        )
    return pd.DataFrame(rows)


def prepare_data(project_root: Path, output_dir: Path, config: Optional[ModelConfig] = None) -> Dict[str, Any]:
    """Create aggregate priors; sequence and strain-specific records are excluded."""

    config = config or ModelConfig()
    output_dir.mkdir(parents=True, exist_ok=True)
    sabio = summarize_sabio(project_root / "data" / "sabio-rk" / "raw")
    sabio.to_csv(output_dir / "sabio_kinetic_priors.csv", index=False)

    brenda_files = list((project_root / "data" / "brenda").glob("*.json"))
    brenda = summarize_brenda(brenda_files[0]) if brenda_files else pd.DataFrame()
    brenda.to_csv(output_dir / "brenda_kinetic_priors.csv", index=False)

    registry = parameter_registry(config)
    registry.to_csv(output_dir / "parameter_registry.csv", index=False)
    summary = {
        "scope": "aggregate physicochemical and population-scale priors only",
        "excluded": ["protein sequences", "strain-specific records", "genetic circuits", "mutation sites"],
        "target_reaction_classes": list(TARGET_EC),
        "sabio_summary_rows": int(len(sabio)),
        "brenda_summary_rows": int(len(brenda)),
        "parameter_registry_rows": int(len(registry)),
    }
    (output_dir / "data_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary

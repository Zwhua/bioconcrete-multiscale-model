"""Carbonate speciation and a reproducible geochemical lookup table."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import shutil
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

from .config import ChemistryConfig


def carbonate_fractions(ph: np.ndarray, temperature_c: float = 25.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return CO2(aq), bicarbonate, and carbonate fractions of total inorganic C.

    The temperature correction is intentionally modest. The lookup metadata records
    when this analytical fallback is used instead of a PHREEQC executable.
    """

    ph_values = np.asarray(ph, dtype=float)
    delta = temperature_c - 25.0
    pka1 = 6.35 - 0.005 * delta
    pka2 = 10.33 - 0.010 * delta
    hydrogen = np.power(10.0, -ph_values)
    ka1 = 10.0 ** (-pka1)
    ka2 = 10.0 ** (-pka2)
    denominator = hydrogen * hydrogen + ka1 * hydrogen + ka1 * ka2
    return (
        hydrogen * hydrogen / denominator,
        ka1 * hydrogen / denominator,
        ka1 * ka2 / denominator,
    )


def calcite_saturation(
    calcium_mol_m3: np.ndarray,
    inorganic_carbon_mol_m3: np.ndarray,
    ph: np.ndarray,
    chemistry: ChemistryConfig,
    temperature_c: float = 25.0,
) -> np.ndarray:
    """Approximate calcite saturation ratio using molar activities."""

    _, _, alpha_co3 = carbonate_fractions(np.asarray(ph), temperature_c)
    ca_molar = np.maximum(np.asarray(calcium_mol_m3), 0.0) / 1000.0
    co3_molar = np.maximum(np.asarray(inorganic_carbon_mol_m3), 0.0) * alpha_co3 / 1000.0
    ionic_strength_correction = 1.0 / (1.0 + 0.15 * np.sqrt(np.maximum(ca_molar, 0.0)))
    return np.maximum(ca_molar * co3_molar * ionic_strength_correction**2 / chemistry.calcite_ksp, 0.0)


class GeochemLookup:
    """Interpolated geochemical table with an analytical fallback."""

    def __init__(self, frame: Optional[pd.DataFrame] = None):
        self.frame = frame
        self._interpolator = None
        self._bounds = None
        if frame is not None and not frame.empty:
            ph = np.sort(frame["ph"].unique())
            ca = np.sort(frame["calcium_mol_m3"].unique())
            ct = np.sort(frame["inorganic_carbon_mol_m3"].unique())
            ordered = frame.set_index(["ph", "calcium_mol_m3", "inorganic_carbon_mol_m3"]).sort_index()
            values = ordered["calcite_saturation"].to_numpy().reshape(len(ph), len(ca), len(ct))
            self._interpolator = RegularGridInterpolator((ph, ca, ct), values, bounds_error=False, fill_value=None)
            self._bounds = ((ph[0], ph[-1]), (ca[0], ca[-1]), (ct[0], ct[-1]))

    @classmethod
    def load(cls, path: Path) -> "GeochemLookup":
        return cls(pd.read_csv(path))

    def saturation(
        self,
        calcium: np.ndarray,
        carbon: np.ndarray,
        ph: float,
        chemistry: ChemistryConfig,
        temperature_c: float,
    ) -> np.ndarray:
        ca = np.asarray(calcium, dtype=float)
        ct = np.asarray(carbon, dtype=float)
        if self._interpolator is None:
            return calcite_saturation(ca, ct, np.full_like(ca, ph), chemistry, temperature_c)
        points = np.column_stack([np.full(ca.size, ph), ca.ravel(), ct.ravel()])
        if self._bounds is not None:
            for index, (lower, upper) in enumerate(self._bounds):
                points[:, index] = np.clip(points[:, index], lower, upper)
        values = self._interpolator(points).reshape(ca.shape)
        return np.nan_to_num(np.maximum(values, 0.0), nan=0.0, posinf=1e12, neginf=0.0)


def _contains_phase(database: Path, names: Iterable[str]) -> Dict[str, bool]:
    if not database.exists():
        return {name: False for name in names}
    text = database.read_text(encoding="utf-8", errors="ignore").lower()
    return {name: name.lower() in text for name in names}


def build_geochem_grid(
    output_dir: Path,
    chemistry: ChemistryConfig,
    project_root: Path,
    temperature_c: float = 30.0,
) -> Tuple[Path, Path]:
    """Build the portable lookup used by spatial simulations.

    A local PHREEQC executable is detected and reported, but this release uses the
    transparent carbonate-equilibrium calculation for the grid. This keeps results
    reproducible on Windows while retaining database phase validation.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    ph_values = np.linspace(8.0, 13.0, 21)
    calcium_values = np.unique(np.concatenate(([0.0], np.geomspace(0.01, 15000.0, 30))))
    carbon_values = np.unique(np.concatenate(([0.0], np.geomspace(0.01, 30000.0, 32))))
    rows = []
    for ph in ph_values:
        ca_grid, ct_grid = np.meshgrid(calcium_values, carbon_values, indexing="ij")
        omega = calcite_saturation(ca_grid, ct_grid, np.full_like(ca_grid, ph), chemistry, temperature_c)
        alpha = carbonate_fractions(np.array([ph]), temperature_c)
        for i, ca in enumerate(calcium_values):
            for j, ct in enumerate(carbon_values):
                rows.append((ph, ca, ct, alpha[0][0], alpha[1][0], alpha[2][0], omega[i, j]))
    grid_path = output_dir / "carbonate_lookup.csv"
    pd.DataFrame(
        rows,
        columns=[
            "ph",
            "calcium_mol_m3",
            "inorganic_carbon_mol_m3",
            "alpha_co2",
            "alpha_hco3",
            "alpha_co3",
            "calcite_saturation",
        ],
    ).to_csv(grid_path, index=False)

    phreeqc_database = project_root / "data" / "phreeqc" / "raw" / "database" / "phreeqc.dat"
    concrete_database = project_root / "data" / "phreeqc" / "raw" / "database" / "Concrete_PHR.dat"
    cemdata_database = project_root / "data" / "cemdata" / "raw" / "CEMDATA18-31-03-2022-phaseVol.dat"
    metadata = {
        "backend": "analytical carbonate equilibrium",
        "phreeqc_executable": shutil.which("phreeqc"),
        "temperature_c": temperature_c,
        "warning": "The analytical grid is a model prior, not a replacement for pore-solution measurements.",
        "databases": {
            "phreeqc": str(phreeqc_database),
            "concrete": str(concrete_database),
            "cemdata": str(cemdata_database),
        },
        "phase_checks": {
            "phreeqc": _contains_phase(phreeqc_database, ["Calcite", "Portlandite"]),
            "concrete": _contains_phase(concrete_database, ["Portlandite", "Tobermorite", "Jennite"]),
            "cemdata": _contains_phase(cemdata_database, ["Cal", "Arg", "Portlandite", "CSH"]),
        },
    }
    metadata_path = output_dir / "geochem_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return grid_path, metadata_path

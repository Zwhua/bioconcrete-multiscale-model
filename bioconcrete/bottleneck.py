"""Legacy heuristic bottleneck classification retained for V4 reproducibility."""

from __future__ import annotations

from typing import Dict, Mapping


def classify_bottleneck(indicators: Mapping[str, float]) -> Dict[str, object]:
    """Return the dominant limitation and auditable normalized scores."""

    scores = {
        "activity_limited": 1.0 - indicators.get("activity_saturation", 1.0),
        "substrate_limited": 1.0 - indicators.get("substrate_saturation", 1.0),
        "calcium_limited": 1.0 - indicators.get("calcium_saturation", 1.0),
        "oxygen_limited": 1.0 - indicators.get("oxygen_saturation", 1.0),
        "transport_limited": 1.0 - indicators.get("transport_access", 1.0),
        "geometry_limited": 1.0 - indicators.get("geometry_efficiency", 1.0),
        "release_limited": 1.0 - indicators.get("release_fraction", 1.0),
    }
    scores = {name: min(max(float(value), 0.0), 1.0) for name, value in scores.items()}
    dominant = max(scores, key=scores.get)
    return {
        "dominant_bottleneck": dominant, "score": scores[dominant],
        "scores": scores,
        "basis": "legacy_heuristic; excluded from V5 scientific conclusions",
        "evidence_class": "legacy_heuristic",
    }


def indicators_from_row(row: Mapping[str, float]) -> Dict[str, float]:
    """Derive bounded indicators from a simulation/design-matrix row."""

    width = max(float(row.get("crack_width_mm", 0.3)), 1e-9)
    wet = float(row.get("wet_hours_per_day", 12.0))
    activity = float(row.get("activity_multiplier", 1.0))
    dosage = float(row.get("agent_dosage", 1.0))
    return {
        "activity_saturation": min(activity / 2.0, 1.0),
        "substrate_saturation": min(dosage, 1.0),
        "calcium_saturation": min(dosage, 1.0),
        "oxygen_saturation": min(wet / 12.0, 1.0),
        "transport_access": min(wet / 24.0 + 0.25, 1.0),
        "geometry_efficiency": min(0.3 / width, 1.0),
        "release_fraction": min(float(row.get("release_fraction", 0.5)), 1.0),
    }

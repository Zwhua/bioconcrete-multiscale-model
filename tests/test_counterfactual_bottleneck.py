import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from bioconcrete.counterfactual import (
    control_coefficient,
    counterfactual_bottleneck,
    summarize_controls,
)


class CounterfactualBottleneckTests(unittest.TestCase):
    def test_control_coefficient_uses_model_response(self):
        self.assertAlmostEqual(control_coefficient(2.0, 2.4, 0.1), 2.0)

    def test_parallel_factors_are_reported(self):
        frame = pd.DataFrame([
            {"factor": factor, "perturbation_fraction": delta,
             "control_closure_28d": value, "control_calcite_mass_mg": value,
             "control_transmissivity_ratio": value}
            for factor, value in (("activity", 1.0), ("release", 0.95))
            for delta in (-0.1, 0.1)
        ])
        _, summary = summarize_controls(frame)
        self.assertTrue(summary["dominant_bottleneck"].startswith("parallel:"))

    def test_workflow_writes_required_outputs(self):
        calls = {"count": 0}

        def fake_output(config):
            calls["count"] += 1
            return {"closure_28d": config.kinetics.activity_multiplier,
                    "calcite_mass_mg": config.kinetics.capsule_release_s,
                    "transmissivity_ratio": config.transport.crack_width_mm}

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch("bioconcrete.counterfactual._outputs", side_effect=fake_output):
                result = counterfactual_bottleneck(output, perturbations=(-0.1, 0.1))
            self.assertEqual(result["basis"], "normalized counterfactual model responses")
            for name in (
                "counterfactual_control_coefficients.csv", "dominant_bottlenecks.csv",
                "bottleneck_uncertainty.csv", "bottleneck_summary.json",
            ):
                self.assertTrue((output / name).exists())
            self.assertEqual(calls["count"], 1 + 2 * 8)


if __name__ == "__main__":
    unittest.main()

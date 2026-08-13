import json
import tempfile
import unittest
from pathlib import Path

from bioconcrete.config import ModelConfig
from bioconcrete.dimensionless import analyze_dimensionless, write_dimensionless_summary


class DimensionlessTests(unittest.TestCase):
    def test_definitions_and_default_screen(self):
        config = ModelConfig()
        result = analyze_dimensionless(config)
        expected_y = (config.transport.crack_width_mm * 1e-3) ** 2 / config.transport.diffusivity_oxygen_m2_s
        self.assertAlmostEqual(result.diffusion_time_aperture_s, expected_y)
        self.assertAlmostEqual(result.peclet_length, 0.0)
        self.assertIsNone(result.advection_time_length_s)
        self.assertAlmostEqual(result.damkohler_length, config.kinetics.maximum_growth_s * (0.1 ** 2) / config.transport.diffusivity_oxygen_m2_s)
        self.assertTrue(result.two_point_five_d_applicable)

    def test_slow_aperture_mixing_requires_3d(self):
        config = ModelConfig()
        result = analyze_dimensionless(config, diffusivity_m2_s=1e-14)
        self.assertFalse(result.two_point_five_d_applicable)
        self.assertIn("full 3D", result.applicability_reason)

    def test_summary_is_atomically_serializable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "dimensionless.json"
            write_dimensionless_summary(analyze_dimensionless(ModelConfig()), path)
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("Infinity", text)
            payload = json.loads(text)
            self.assertEqual(payload["evidence_label"], "uncalibrated 3D model output; not experimental data")
            self.assertFalse(path.with_suffix(".json.tmp").exists())


if __name__ == "__main__":
    unittest.main()

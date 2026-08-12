import unittest

from bioconcrete.bottleneck import classify_bottleneck


class BottleneckTests(unittest.TestCase):
    def test_geometry_limitation_is_detected(self):
        result = classify_bottleneck({
            "activity_saturation": 1, "substrate_saturation": 1, "calcium_saturation": 1,
            "oxygen_saturation": 1, "transport_access": 1, "geometry_efficiency": 0.1,
            "release_fraction": 1,
        })
        self.assertEqual(result["dominant_bottleneck"], "geometry_limited")
        self.assertNotIn("experiment", result["evidence_class"])


if __name__ == "__main__":
    unittest.main()

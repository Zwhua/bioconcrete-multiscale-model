import json
import unittest
from pathlib import Path


class LegacyRegressionFixtureTests(unittest.TestCase):
    def test_frozen_v051_artifacts_match_fixture(self):
        fixture=json.loads((Path('tests/fixtures/v0_5_1_regression.json')).read_text(encoding='utf-8'))
        paths={'0d':'20260811_231012_0d','1d':'20260811_231144_1d','2d':'20260811_231655_2d'}
        for level,directory in paths.items():
            actual=json.loads((Path('model_runs/final')/directory/'summary.json').read_text(encoding='utf-8'))
            for key,expected in fixture['levels'][level].items():
                self.assertAlmostEqual(actual[key],expected,delta=abs(expected)*fixture['rtol'])
        self.assertEqual(fixture['source_commit'][:7],'1890526')

if __name__ == '__main__': unittest.main()

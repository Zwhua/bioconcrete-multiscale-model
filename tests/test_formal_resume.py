import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from bioconcrete.config import ModelConfig
from bioconcrete.formal_analysis import _sample_id, evaluate_resumable


class FormalResumeTests(unittest.TestCase):
    def test_resume_does_not_repeat_samples(self):
        problem = {"names": ["kinetics.activity_multiplier"]}
        matrices = {"morris": np.array([[1.0], [2.0]]), "sobol": np.array([[0.5], [1.5]])}
        calls = {"count": 0}

        def fake(task):
            calls["count"] += 1
            method, index, row, _, _ = task
            return {"sample_id": _sample_id(method, index, row), "method": method,
                    "sample_index": index, "kinetics.activity_multiplier": row[0],
                    "response": row[0]}

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch("bioconcrete.formal_analysis._evaluate_sample", side_effect=fake):
                first, _ = evaluate_resumable(matrices, problem, ModelConfig(), output)
                initial_calls = calls["count"]
                second, _ = evaluate_resumable(matrices, problem, ModelConfig(), output, resume=True)
            self.assertEqual(initial_calls, 4)
            self.assertEqual(calls["count"], initial_calls)
            np.testing.assert_allclose(first["sobol"], second["sobol"])


if __name__ == "__main__":
    unittest.main()

import unittest

import numpy as np

from bioconcrete.config import ModelConfig
from bioconcrete.model import _initial_state, _reaction_step
from bioconcrete.reaction_kernel import reaction_step_cells
from bioconcrete.state import STATE_NAMES


class ReactionKernelTests(unittest.TestCase):
    def setUp(self):
        self.config = ModelConfig()
        self.state = _initial_state(self.config, np.ones(5), "0d")

    def test_model_reexports_canonical_schema(self):
        from bioconcrete import model
        self.assertIs(model.STATE_NAMES, STATE_NAMES)

    def test_batches_match_existing_single_source(self):
        expected = _reaction_step(self.state, 0.0, 60.0, self.config, None)
        observed = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=2)
        np.testing.assert_allclose(observed, expected, rtol=2e-5, atol=1e-9)

    def test_batch_size_changes_only_within_solver_tolerance(self):
        first = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=1)
        second = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=4)
        np.testing.assert_allclose(first, second, rtol=2e-5, atol=1e-9)

    def test_invalid_shape_is_rejected(self):
        with self.assertRaises(ValueError):
            reaction_step_cells(np.ones((2, 3)), 0.0, 1.0, self.config)

    def test_serial_and_parallel_are_identical(self):
        serial = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=2, workers=1)
        parallel = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=2, workers=2,
                                       parallel_backend="thread")
        np.testing.assert_array_equal(serial, parallel)

    def test_serial_and_process_are_identical(self):
        serial = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=2)
        process = reaction_step_cells(self.state, 0.0, 60.0, self.config, batch_size=2, workers=2,
                                      parallel_backend="process")
        np.testing.assert_array_equal(serial, process)


if __name__ == "__main__":
    unittest.main()

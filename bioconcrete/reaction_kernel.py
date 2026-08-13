"""Shared local reaction integration, including deterministic cell batching.

The scientific equations remain implemented once in :mod:`bioconcrete.model`
during the compatibility refactor. This module is the dimension-independent
entry point used by new solvers and deliberately delegates to that source.
"""

from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import List, Optional, Tuple

import numpy as np

from .chemistry import GeochemLookup
from .config import ModelConfig
from .state import S, STATE_NAMES


@dataclass(frozen=True)
class ReactionBatchDiagnostics:
    cell_count: int
    batch_size: int
    batch_count: int
    failed_cell_indices: List[int]


def reaction_step_cells(
    state: np.ndarray,
    time_s: float,
    dt_s: float,
    config: ModelConfig,
    geochem: Optional[GeochemLookup] = None,
    *,
    batch_size: int = 64,
    workers: int = 1,
    parallel_backend: str = "serial",
) -> np.ndarray:
    """Integrate independent reaction cells in deterministic contiguous batches."""
    from .model import _reaction_step

    values = np.asarray(state, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(STATE_NAMES):
        raise ValueError("state must have shape (n_cells, {})".format(len(STATE_NAMES)))
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if workers <= 0:
        raise ValueError("workers must be positive")
    if parallel_backend not in {"serial", "thread", "process"}:
        raise ValueError("parallel_backend must be serial, thread, or process")
    updated = np.empty_like(values)
    ranges = [(start, min(start + batch_size, values.shape[0]))
              for start in range(0, values.shape[0], batch_size)]
    tasks = [(start, stop, values[start:stop], time_s, dt_s, config, geochem)
             for start, stop in ranges]
    executor = None
    if workers == 1 or len(ranges) <= 1 or parallel_backend == "serial":
        results = map(_integrate_reaction_batch, tasks)
    elif parallel_backend == "thread":
        executor = ThreadPoolExecutor(max_workers=workers)
        results = executor.map(_integrate_reaction_batch, tasks)
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        results = executor.map(_integrate_reaction_batch, tasks, chunksize=1)
    try:
        for start, stop, result in results:
            updated[start:stop] = result
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
    return updated


def _integrate_reaction_batch(task):
    """Pickle-safe worker used by thread and process executors."""
    from .model import _reaction_step
    start, stop, values, time_s, dt_s, config, geochem = task
    try:
        result = _reaction_step(values, time_s, dt_s, config, geochem)
    except Exception as error:
        indices = list(range(start, stop))
        raise RuntimeError(
            "Reaction batch failed for cell indices {}: {}".format(indices, error)
        ) from error
    return start, stop, result

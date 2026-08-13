"""Multiscale physicochemical model for self-healing concrete."""

from .config import ModelConfig
from .model import SimulationResult, simulate_0d, simulate_1d, simulate_2d

__all__ = ["ModelConfig", "SimulationResult", "simulate_0d", "simulate_1d", "simulate_2d"]
__version__ = "0.5.1"

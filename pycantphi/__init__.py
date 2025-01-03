"""
PyCantPhi: A Python library for calculating anthropogenic carbon in the North Atlantic Ocean.

This package implements the φ CT° method for calculating anthropogenic carbon,
including support for water mass analysis and CO2 system calculations.
"""

from importlib.metadata import version, PackageNotFoundError

# Version and metadata
try:
    __version__ = version("pycantphi")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = "Raphaël Bajon"
__email__ = "raphael.bajon@ifremer.fr"

# Import main classes and functions for convenient access
from .core.calculator import CantPhiCt0 as cantphi
from .core.cant import CantCalculator 

# Define public API
__all__ = [
    'cantphi',
    'CantCalculator'
]

# Optional: Package initialization
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())  # Let application configure logging

# Optional: Default configuration

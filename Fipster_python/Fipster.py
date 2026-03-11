"""
Legacy compatibility shim for projects that still import `Fipster`.
"""

from fipster import FIP_signal, Sweepset, __version__
from fipster.cli import main

__all__ = ["FIP_signal", "Sweepset", "__version__", "main"]

"""Processors: the abstract :class:`Processor` base class and concrete ones.

Concrete processors (added over time) live in their own modules in this
subpackage and subclass :class:`Processor`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mf_engine.processors.base import Processor

__all__ = [
    "Processor",
    "GNPSProcessor",
    "GNPSInput",
    "GNPSResult",
    "ANNOTATION_COLUMNS",
]

# The GNPS symbols exposed lazily by __getattr__ (all live in the gnps module).
# ANNOTATION_COLUMNS is the public output schema — the ordered list of columns in
# every GNPSResult record — exposed so callers can build headers / validate.
_LAZY_GNPS_SYMBOLS = ("GNPSProcessor", "GNPSInput", "GNPSResult", "ANNOTATION_COLUMNS")

if TYPE_CHECKING:  # for type checkers / IDEs only
    from mf_engine.processors.gnps import (
        ANNOTATION_COLUMNS,
        GNPSInput,
        GNPSProcessor,
        GNPSResult,
    )


def __getattr__(name: str) -> Any:
    """Lazily expose the pandas-backed GNPS symbols (PEP 562).

    This keeps ``import mf_engine.processors`` dependency-free; pandas is only
    imported when a GNPS symbol is actually accessed.
    """
    if name in _LAZY_GNPS_SYMBOLS:
        from mf_engine.processors import gnps

        return getattr(gnps, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

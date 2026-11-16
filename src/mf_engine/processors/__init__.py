"""Processors: the abstract :class:`Processor` base class and concrete ones.

Concrete processors (added over time) live in their own modules in this
subpackage and subclass :class:`Processor`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from mf_engine.processors.base import Processor

__all__ = ["Processor", "GNPSProcessor", "GNPSInput", "GNPSResult"]

if TYPE_CHECKING:  # for type checkers / IDEs only
    from mf_engine.processors.gnps import GNPSInput, GNPSProcessor, GNPSResult


def __getattr__(name: str) -> Any:
    """Lazily expose the pandas-backed GNPS symbols (PEP 562).

    This keeps ``import mf_engine.processors`` dependency-free; pandas is only
    imported when a GNPS symbol is actually accessed.
    """
    if name in ("GNPSProcessor", "GNPSInput", "GNPSResult"):
        from mf_engine.processors import gnps

        return getattr(gnps, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

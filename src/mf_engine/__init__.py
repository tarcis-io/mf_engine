"""mf_engine — the MassFlow engine.

A reusable, importable library of *processors* for MZmine + GNPS metabolomics
batches. The public surface is intentionally small:

    from mf_engine import Processor

Concrete processors subclass :class:`Processor` and implement ``process``.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

# Only the dependency-free base class is imported eagerly. Concrete processors
# that need third-party libraries (e.g. GNPSProcessor needs pandas) are imported
# from their own modules, so ``import mf_engine`` stays lightweight:
#     from mf_engine.processors.gnps import GNPSProcessor
from mf_engine.processors.base import Processor

__all__ = ["Processor", "__version__"]

# Single source of truth: the version declared in pyproject.toml, read from the
# installed package metadata. Falls back to "0.0.0" when running from a source
# tree that has not been installed.
try:
    __version__ = version("mf_engine")
except PackageNotFoundError:  # pragma: no cover - source tree, not installed
    __version__ = "0.0.0"

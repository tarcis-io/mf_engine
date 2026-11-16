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

# Version reported when the package metadata is unavailable (i.e. running from a
# source tree that was never installed).
_FALLBACK_VERSION = "0.0.0"


def _read_version() -> str:
    """Return the installed package version, or the source-tree fallback.

    Single source of truth is the version declared in ``pyproject.toml``, read
    from the installed distribution metadata. When the package is imported from a
    source tree that has not been installed, the metadata is unavailable and this
    returns :data:`_FALLBACK_VERSION` rather than raising.
    """
    try:
        return version("mf_engine")
    except PackageNotFoundError:
        return _FALLBACK_VERSION


__version__ = _read_version()

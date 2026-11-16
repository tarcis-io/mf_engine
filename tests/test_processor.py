"""Tests for the Processor abstract base class and package surface."""

from __future__ import annotations

import importlib

import pytest

import mf_engine
from mf_engine import Processor


def test_package_exposes_version():
    # A string version is always available (from metadata, or the source-tree
    # fallback), so callers can rely on mf_engine.__version__.
    assert isinstance(mf_engine.__version__, str)
    assert mf_engine.__version__


def test_version_falls_back_when_not_installed(monkeypatch):
    # When the distribution metadata is unavailable (running from a source tree
    # that was never installed), __version__ falls back to "0.0.0" rather than
    # raising on import.
    from importlib.metadata import PackageNotFoundError

    def _raise(_name: str) -> str:
        raise PackageNotFoundError

    monkeypatch.setattr("importlib.metadata.version", _raise)
    reloaded = importlib.reload(mf_engine)
    try:
        assert reloaded.__version__ == "0.0.0"
    finally:
        # Restore real metadata-backed version for any later tests.
        monkeypatch.undo()
        importlib.reload(mf_engine)


def test_processor_is_abstract():
    # The ABC cannot be instantiated directly.
    with pytest.raises(TypeError):
        Processor()  # type: ignore[abstract]


def test_subclass_must_implement_process():
    class Incomplete(Processor):
        pass

    with pytest.raises(TypeError):
        Incomplete()  # type: ignore[abstract]


def test_concrete_subclass_works():
    class Upper(Processor):
        def process(self, data: str) -> str:
            return data.upper()

    proc = Upper()
    assert proc.process("hello") == "HELLO"


def test_default_name_is_class_name():
    class Foo(Processor):
        def process(self):  # pragma: no cover - trivial
            return None

    assert Foo().name == "Foo"
    assert Foo(name="custom").name == "custom"


def test_repr():
    class Bar(Processor):
        def process(self):  # pragma: no cover - trivial
            return None

    assert repr(Bar(name="x")) == "Bar(name='x')"


def test_processors_lazily_exposes_gnps_symbols():
    # The processors subpackage stays import-light (no pandas at import time) and
    # exposes the GNPS symbols on demand via PEP 562 __getattr__.
    from mf_engine import processors

    for name in ("GNPSProcessor", "GNPSInput", "GNPSResult"):
        assert getattr(processors, name) is not None


def test_processors_unknown_attr_raises_attribute_error():
    # The lazy __getattr__ must not swallow genuine typos / missing names.
    from mf_engine import processors

    with pytest.raises(AttributeError):
        processors.NoSuchSymbol  # noqa: B018 - intentional attribute access

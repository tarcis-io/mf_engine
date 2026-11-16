"""The abstract base class for all MassFlow processors.

A *processor* takes some input and produces some output. The base class fixes
only the minimal contract shared by every processor — an optional name and a
single abstract :meth:`process` method — leaving the concrete input/output
types to each subclass.

Example
-------
    from mf_engine import Processor

    class UppercaseProcessor(Processor):
        def process(self, data: str) -> str:
            return data.upper()

    proc = UppercaseProcessor()
    proc.process("hello")  # -> "HELLO"
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class Processor(ABC):
    """Abstract base class for a processing unit.

    Subclasses must implement :meth:`process`. The base class provides a shared
    initialiser that records an optional ``name`` (defaulting to the subclass
    name), which is handy for logging and identifying processors at run time.
    """

    def __init__(self, name: str | None = None) -> None:
        self.name: str = name or type(self).__name__

    @abstractmethod
    def process(self, *args: Any, **kwargs: Any) -> Any:
        """Run the processor's work and return its result.

        Concrete subclasses define the exact signature and return type.
        """
        raise NotImplementedError

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r})"

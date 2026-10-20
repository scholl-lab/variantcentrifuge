"""Accumulate wall-clock seconds per named association sub-step."""

from __future__ import annotations

import logging
import time
from collections.abc import Iterator
from contextlib import contextmanager

logger = logging.getLogger(__name__)


@contextmanager
def substep(timings: dict[str, float], name: str) -> Iterator[None]:
    """Add the wall-clock seconds of the ``with`` body to ``timings[name]``."""
    start = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        timings[name] = timings.get(name, 0.0) + elapsed
        logger.debug("Association sub-step %s: %.3f s", name, elapsed)

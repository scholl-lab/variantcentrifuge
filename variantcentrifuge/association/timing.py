"""Accumulate wall-clock seconds per named association sub-step."""

from __future__ import annotations

import logging
import resource
import time
from collections.abc import Callable, Iterator
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


class Progress:
    """Throttled progress line (done/total, rate, ETA), logged by the parent process only."""

    def __init__(
        self,
        label: str,
        total: int,
        interval: float = 30.0,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self.label, self.total, self.interval, self._clock = label, total, interval, clock
        self._start = self._last = clock()

    def update(self, done: int) -> None:
        """Log when ``interval`` seconds passed since the last line, and at completion."""
        now = self._clock()
        if done >= self.total or now - self._last >= self.interval:
            self._last = now
            elapsed = max(now - self._start, 1e-9)
            rate = done / elapsed
            eta = int((self.total - done) / rate) if rate else 0
            logger.info(
                "%s: %d/%d (%.1f%%), %.2f/s, ETA %d:%02d:%02d",
                self.label, done, self.total, 100.0 * done / max(self.total, 1), rate,
                eta // 3600, eta % 3600 // 60, eta % 60,
            )  # fmt: skip


def log_stage_summary(timings: dict[str, float], facts: list[str]) -> None:
    """Log one INFO line per sub-step, the total, peak RSS (self / largest child) and ``facts``."""
    for name, seconds in timings.items():
        logger.info("Association sub-step %s: %.1f s", name, seconds)
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024  # KiB on Linux -> MiB
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    logger.info(
        "Association total: %.1f s, peak RSS %.0f MiB (largest worker %.0f MiB)",
        sum(timings.values()), own, child,
    )  # fmt: skip
    for fact in facts:
        logger.info("Association: %s", fact)

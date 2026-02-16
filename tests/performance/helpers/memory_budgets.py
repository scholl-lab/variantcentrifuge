"""
Memory tracking and budget enforcement for performance benchmarks.

Uses tracemalloc to measure peak memory usage and provides warning-only
budget enforcement (never raises exceptions).
"""

import tracemalloc
import warnings

# Default memory budgets (MB) at 2x expected peak for 10K variants x 100 samples
# These will be refined once actual profiling data exists
INHERITANCE_BUDGET_MB = 512
COMP_HET_BUDGET_MB = 256
GENOTYPE_REPLACEMENT_BUDGET_MB = 1024
GENE_BURDEN_BUDGET_MB = 512
SCORING_BUDGET_MB = 256


class MemoryTracker:
    """
    Context manager for tracking memory usage with tracemalloc.

    Captures peak memory usage during execution and provides
    MB-denominated access to current and peak values.

    Examples
    --------
    >>> with MemoryTracker() as tracker:
    ...     data = [i for i in range(1000000)]
    >>> tracker.peak_mb > 0
    True
    """

    def __init__(self):
        """Initialize the memory tracker."""
        self._peak = 0
        self._current = 0

    def __enter__(self):
        """Start tracemalloc and reset peak."""
        tracemalloc.start()
        tracemalloc.reset_peak()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Capture peak memory and stop tracemalloc."""
        current, peak = tracemalloc.get_traced_memory()
        self._current = current
        self._peak = peak
        tracemalloc.stop()
        return False

    @property
    def peak_mb(self) -> float:
        """Peak memory usage in megabytes."""
        return self._peak / (1024 * 1024)

    @property
    def current_mb(self) -> float:
        """Current memory usage in megabytes."""
        return self._current / (1024 * 1024)


def warn_if_over_budget(peak_mb: float, budget_mb: float, context: str) -> None:
    """
    Issue warning if peak memory exceeds budget.

    Never raises exceptions - memory violations are warning-only per
    project guidelines. This allows benchmarks to complete while
    flagging concerning memory usage for investigation.

    Parameters
    ----------
    peak_mb : float
        Peak memory usage in megabytes
    budget_mb : float
        Budget threshold in megabytes
    context : str
        Description of what was being benchmarked (for warning message)

    Examples
    --------
    >>> warn_if_over_budget(100.0, 50.0, "test operation")
    # Issues UserWarning
    >>> warn_if_over_budget(50.0, 100.0, "test operation")
    # No warning
    """
    if peak_mb > budget_mb:
        warnings.warn(
            f"Memory budget exceeded for {context}: "
            f"peak={peak_mb:.1f} MB, budget={budget_mb} MB "
            f"(over by {peak_mb - budget_mb:.1f} MB)",
            UserWarning,
            stacklevel=2,
        )

"""
Pytest fixtures for performance benchmarks.

Provides factory fixtures for synthetic data generation and memory tracking.
All data is fully synthetic with reproducible seeding.
"""

import pytest

from .helpers.memory_budgets import MemoryTracker
from .helpers.synthetic_data import (
    generate_synthetic_gene_burden_data,
    generate_synthetic_pedigree,
    generate_synthetic_scoring_config,
    generate_synthetic_variants,
)


@pytest.fixture
def synthetic_variants():
    """
    Factory fixture for generating synthetic variant DataFrames.

    Returns a callable that generates DataFrames matching variantcentrifuge
    structure with configurable size and reproducible seeding.

    Returns
    -------
    callable
        Function(n_variants, n_samples, seed=42) -> pd.DataFrame

    Examples
    --------
    >>> def test_something(synthetic_variants):
    ...     df = synthetic_variants(100, 10)
    ...     assert df.shape == (100, 9)
    ...     # Test with different seed
    ...     df2 = synthetic_variants(100, 10, seed=123)
    ...     assert not df.equals(df2)
    ...     # Same seed gives same data
    ...     df3 = synthetic_variants(100, 10, seed=123)
    ...     assert df2.equals(df3)
    """

    def _generate(n_variants: int, n_samples: int, seed: int = 42):
        return generate_synthetic_variants(n_variants, n_samples, seed)

    return _generate


@pytest.fixture
def synthetic_pedigree():
    """
    Factory fixture for generating synthetic pedigree data.

    Returns a callable that generates pedigree dictionaries matching
    variantcentrifuge format with configurable size and reproducible seeding.

    Returns
    -------
    callable
        Function(n_samples=3, seed=42) -> dict

    Examples
    --------
    >>> def test_something(synthetic_pedigree):
    ...     ped = synthetic_pedigree(10)
    ...     assert len(ped) == 10
    ...     assert 'CHILD_001' in ped
    """

    def _generate(n_samples: int = 3, seed: int = 42):
        return generate_synthetic_pedigree(n_samples, seed)

    return _generate


@pytest.fixture
def synthetic_scoring_config():
    """
    Fixture for generating synthetic scoring configuration.

    Returns minimal scoring config matching variantcentrifuge format.
    Config is deterministic (not a factory - always the same).

    Returns
    -------
    dict
        Scoring configuration with variables and formulas

    Examples
    --------
    >>> def test_something(synthetic_scoring_config):
    ...     assert 'variables' in synthetic_scoring_config
    ...     assert 'formulas' in synthetic_scoring_config
    """
    return generate_synthetic_scoring_config()


@pytest.fixture
def synthetic_gene_burden_data():
    """
    Factory fixture for generating synthetic gene burden test data.

    Returns a callable that generates DataFrames with gene burden GT format
    plus case/control sample sets.

    Returns
    -------
    callable
        Function(n_variants, n_samples, n_genes=50, seed=42) -> tuple[pd.DataFrame, set, set]

    Examples
    --------
    >>> def test_something(synthetic_gene_burden_data):
    ...     df, cases, controls = synthetic_gene_burden_data(100, 20, n_genes=10)
    ...     assert df.shape[0] == 100
    ...     assert len(cases) + len(controls) == 20
    """

    def _generate(n_variants: int, n_samples: int, n_genes: int = 50, seed: int = 42):
        return generate_synthetic_gene_burden_data(n_variants, n_samples, n_genes, seed)

    return _generate


@pytest.fixture
def memory_tracker():
    """
    Fixture for memory tracking with tracemalloc.

    Returns a fresh MemoryTracker instance that can be used as a context manager.

    Returns
    -------
    MemoryTracker
        Context manager for tracking peak memory usage

    Examples
    --------
    >>> def test_something(memory_tracker):
    ...     with memory_tracker as tracker:
    ...         data = [i for i in range(1000000)]
    ...     assert tracker.peak_mb > 0
    """
    return MemoryTracker()

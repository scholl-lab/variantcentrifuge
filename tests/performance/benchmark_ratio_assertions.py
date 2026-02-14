"""
Ratio assertion tests comparing vectorized vs sequential implementations.

These tests measure both implementations in the same run and compute speedup
ratios. Zero flakiness because both implementations run on the same machine
in the same test under identical conditions.

Uses time.perf_counter() directly (NOT pytest-benchmark) to enable within-test
timing comparisons.
"""

import time

import pytest

from variantcentrifuge.inheritance.comp_het import analyze_gene_for_compound_het
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)


def _filter_to_single_gene(df, min_variants=5):
    """
    Filter DataFrame to a single gene with multiple variants.

    Parameters
    ----------
    df : pd.DataFrame
        Variant DataFrame with GENE column
    min_variants : int
        Minimum number of variants required in the gene

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only variants from one gene
    """
    # Find gene with most variants
    gene_counts = df["GENE"].value_counts()
    gene_with_most = gene_counts.index[0]
    gene_df = df[df["GENE"] == gene_with_most].copy()

    # Verify sufficient variants
    if len(gene_df) < min_variants:
        pytest.skip(
            f"Generated data has insufficient variants in largest gene "
            f"({len(gene_df)} < {min_variants})"
        )

    return gene_df


def _time_function_iterations(func, *args, iterations=5, **kwargs):
    """
    Time a function over multiple iterations and return mean time.

    Parameters
    ----------
    func : callable
        Function to time
    args : tuple
        Positional arguments to func
    iterations : int
        Number of iterations to run
    kwargs : dict
        Keyword arguments to func

    Returns
    -------
    float
        Mean execution time in seconds
    """
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        func(*args, **kwargs)
        end = time.perf_counter()
        times.append(end - start)

    return sum(times) / len(times)


@pytest.mark.performance
def test_comp_het_vectorized_vs_original_ratio(synthetic_variants, synthetic_pedigree):
    """
    Compare vectorized vs sequential compound het on a single gene.

    Generates a dataset, filters to one gene, times both implementations,
    and asserts that vectorized is faster (speedup > 1.0).
    """
    # Generate data: 1000 variants, trio pedigree
    df = synthetic_variants(n_variants=1000, n_samples=3, seed=42)
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())

    # Filter to single gene with multiple variants
    gene_df = _filter_to_single_gene(df, min_variants=5)

    print(f"\nTesting compound het ratio on gene with {len(gene_df)} variants")

    # Time sequential implementation
    sequential_mean = _time_function_iterations(
        analyze_gene_for_compound_het,
        gene_df,
        pedigree,
        sample_list,
        iterations=5,
    )
    print(f"Sequential (original): {sequential_mean * 1000:.2f} ms (mean of 5 iterations)")

    # Time vectorized implementation
    vectorized_mean = _time_function_iterations(
        analyze_gene_for_compound_het_vectorized,
        gene_df,
        pedigree,
        sample_list,
        iterations=5,
    )
    print(f"Vectorized: {vectorized_mean * 1000:.2f} ms (mean of 5 iterations)")

    # Compute speedup
    speedup = sequential_mean / vectorized_mean
    print(f"Speedup: {speedup:.2f}x")

    # Assert vectorized is faster
    assert speedup > 1.0, f"Vectorized comp_het should be faster, got {speedup:.2f}x"


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [500, 2000])
def test_comp_het_ratio_at_scale(synthetic_variants, synthetic_pedigree, n_variants):
    """
    Compare vectorized vs sequential compound het at different scales.

    Parametrized to show how speedup ratio changes with dataset size.
    Vectorization benefits should increase with larger datasets.
    """
    # Generate data
    df = synthetic_variants(n_variants=n_variants, n_samples=3, seed=42)
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())

    # Filter to single gene
    gene_df = _filter_to_single_gene(df, min_variants=5)

    print(
        f"\nTesting compound het ratio on {n_variants} variants (gene has {len(gene_df)} variants)"
    )

    # Time both implementations
    sequential_mean = _time_function_iterations(
        analyze_gene_for_compound_het,
        gene_df,
        pedigree,
        sample_list,
        iterations=5,
    )

    vectorized_mean = _time_function_iterations(
        analyze_gene_for_compound_het_vectorized,
        gene_df,
        pedigree,
        sample_list,
        iterations=5,
    )

    # Compute and display speedup
    speedup = sequential_mean / vectorized_mean
    print(
        f"Sequential: {sequential_mean * 1000:.2f} ms, Vectorized: {vectorized_mean * 1000:.2f} ms"
    )
    print(f"Speedup: {speedup:.2f}x")

    # Assert vectorized is faster
    # At small scales (<100 variants), vectorization overhead may dominate
    # We parametrize with n_variants >= 500 where vectorization should win
    assert speedup > 1.0, (
        f"Vectorized comp_het should be faster at {n_variants} variants, "
        f"got {speedup:.2f}x. If this fails, vectorization may not provide "
        f"expected benefits at this scale."
    )

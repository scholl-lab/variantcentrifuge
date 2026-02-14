"""
Inheritance analysis micro and meso benchmarks.

Tests inheritance pattern deduction at single-variant (micro) and full-dataset (meso)
granularity levels. Parametrized across 100, 1K, 10K variant scales.
"""

import pandas as pd
import pytest

from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.inheritance.deducer import deduce_patterns_for_variant


def _expand_gt_column_to_samples(df: pd.DataFrame, sample_list: list[str]) -> pd.DataFrame:
    """
    Transform comma-separated GT column into per-sample columns.

    The inheritance analyzer expects variant_row dict with per-sample columns:
    {"CHILD_001": "0/1", "FATHER_001": "0/0", ...}

    The synthetic_variants fixture generates:
    GT column with comma-separated values: "0/1,0/0,1/1,..."

    This helper expands GT into individual sample columns.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with GT column containing comma-separated genotypes
    sample_list : list[str]
        List of sample IDs in order matching GT column

    Returns
    -------
    pd.DataFrame
        DataFrame with per-sample columns added, GT column preserved
    """
    df = df.copy()

    # Split GT column into individual sample genotypes
    gt_split = df["GT"].str.split(",", expand=True)

    # Add columns for each sample
    for i, sample_id in enumerate(sample_list):
        if i < gt_split.shape[1]:
            df[sample_id] = gt_split[i]

    return df


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_deduce_patterns_scaling(benchmark, synthetic_variants, synthetic_pedigree, n_variants):
    """
    Benchmark full inheritance analysis (meso-level) at multiple scales.

    Tests the complete analyze_inheritance orchestrator which performs:
    1. Per-variant pattern deduction (via df.apply)
    2. Compound heterozygous analysis (per gene)
    3. Pattern prioritization

    This is a meso-level benchmark of the full inheritance module.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate synthetic data with GT column (comma-separated)
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Transform to per-sample columns format for inheritance analysis
    df = _expand_gt_column_to_samples(df, sample_list)

    # Benchmark the full inheritance analysis
    def run_analysis():
        return analyze_inheritance(df, pedigree, sample_list, use_vectorized_comp_het=True)

    if n_variants >= 10000:
        # Use pedantic mode for expensive benchmarks to control iterations
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["component"] = "inheritance_analysis"

    # Verify result structure
    assert isinstance(result, pd.DataFrame)
    assert "Inheritance_Pattern" in result.columns
    assert len(result) == n_variants


@pytest.mark.performance
def test_deduce_single_variant_micro(benchmark, synthetic_variants, synthetic_pedigree):
    """
    Benchmark single-variant pattern deduction (micro-level).

    Tests deduce_patterns_for_variant function in isolation, which is called
    per-row during the meso-level analysis. This micro benchmark isolates
    the per-variant deduction logic from the orchestration overhead.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate small dataset to extract a realistic variant row
    df = synthetic_variants(100, n_samples, seed=42)
    df = _expand_gt_column_to_samples(df, sample_list)

    # Convert first row to dict (how df.apply provides data to the function)
    variant_row = df.iloc[0].to_dict()

    # Benchmark the single-variant deduction function
    result = benchmark(deduce_patterns_for_variant, variant_row, pedigree, sample_list)

    # Store metadata
    benchmark.extra_info["component"] = "deduce_single_variant"
    benchmark.extra_info["n_samples"] = n_samples

    # Verify result
    assert isinstance(result, list)


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_inheritance_vectorized_vs_original(
    benchmark, synthetic_variants, synthetic_pedigree, n_variants
):
    """
    Benchmark inheritance with vectorized compound het enabled.

    This test uses the same code path as test_deduce_patterns_scaling but
    provides a separate benchmark series for tracking vectorized performance
    independently. Useful for comparing vectorized vs original implementations.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)
    df = _expand_gt_column_to_samples(df, sample_list)

    # Benchmark with vectorized compound het explicitly enabled
    def run_analysis():
        return analyze_inheritance(df, pedigree, sample_list, use_vectorized_comp_het=True)

    if n_variants >= 10000:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["component"] = "inheritance_vectorized"
    benchmark.extra_info["use_vectorized_comp_het"] = True

    # Verify
    assert isinstance(result, pd.DataFrame)
    assert len(result) == n_variants

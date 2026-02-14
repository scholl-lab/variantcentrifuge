"""
Compound heterozygous detection benchmarks.

Dedicated benchmarks for compound het analysis, testing both vectorized and
original implementations at multiple scales. Operates at per-gene granularity.
"""

import pandas as pd
import pytest

from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)


def _expand_gt_column_to_samples(df: pd.DataFrame, sample_list: list[str]) -> pd.DataFrame:
    """
    Transform comma-separated GT column into per-sample columns.

    See benchmark_inheritance.py for detailed explanation.
    """
    df = df.copy()
    gt_split = df["GT"].str.split(",", expand=True)
    for i, sample_id in enumerate(sample_list):
        if i < gt_split.shape[1]:
            df[sample_id] = gt_split[i]
    return df


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_comp_het_vectorized_scaling(benchmark, synthetic_variants, synthetic_pedigree, n_variants):
    """
    Benchmark vectorized compound het detection at per-gene level.

    Compound het analysis operates on a single gene's variants at a time.
    This benchmark filters to one gene and measures the vectorized implementation.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)
    df = _expand_gt_column_to_samples(df, sample_list)

    # Filter to single gene (compound het analyzes one gene at a time)
    first_gene = df["GENE"].iloc[0]
    gene_df = df[df["GENE"] == first_gene].copy()
    n_gene_variants = len(gene_df)

    # Benchmark vectorized compound het
    def run_analysis():
        return analyze_gene_for_compound_het_vectorized(gene_df, pedigree, sample_list)

    if n_gene_variants >= 100:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_gene_variants"] = n_gene_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = 1
    benchmark.extra_info["component"] = "comp_het_vectorized"

    # Verify result (returns dict of sample -> list of variant pairs)
    assert isinstance(result, dict)


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000])
def test_comp_het_original_scaling(benchmark, synthetic_variants, synthetic_pedigree, n_variants):
    """
    Benchmark original (non-vectorized) compound het detection.

    Tests the original implementation for comparison. Skip 10K variants
    as the original implementation is significantly slower.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)
    df = _expand_gt_column_to_samples(df, sample_list)

    # Filter to single gene
    first_gene = df["GENE"].iloc[0]
    gene_df = df[df["GENE"] == first_gene].copy()
    n_gene_variants = len(gene_df)

    # Benchmark compound het
    def run_analysis():
        return analyze_gene_for_compound_het_vectorized(gene_df, pedigree, sample_list)

    if n_gene_variants >= 100:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_gene_variants"] = n_gene_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = 1
    benchmark.extra_info["component"] = "comp_het_original"

    # Verify result
    assert isinstance(result, dict)


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_comp_het_multi_gene(benchmark, synthetic_variants, synthetic_pedigree, n_variants):
    """
    Benchmark compound het across multiple genes.

    Simulates real-world usage where compound het is run per-gene across
    many genes. This measures the overhead of gene-by-gene iteration plus
    the vectorized analysis per gene.
    """
    # Generate trio pedigree
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())
    n_samples = len(sample_list)

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)
    df = _expand_gt_column_to_samples(df, sample_list)

    # Count unique genes
    n_genes = df["GENE"].nunique()

    # Benchmark per-gene compound het analysis
    def run_analysis():
        results = {}
        for gene, gene_df in df.groupby("GENE", observed=True):
            gene_results = analyze_gene_for_compound_het_vectorized(gene_df, pedigree, sample_list)
            results[gene] = gene_results
        return results

    if n_variants >= 10000:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["component"] = "comp_het_multi_gene"

    # Verify result
    assert isinstance(result, dict)
    assert len(result) == n_genes

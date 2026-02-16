"""
Gene burden analysis benchmarks.

Tests gene burden analysis performance at multiple scales with realistic
case/control sample distributions.
"""

import pytest

from variantcentrifuge.gene_burden import perform_gene_burden_analysis
from variantcentrifuge.helpers import assign_case_control_counts


@pytest.mark.performance
@pytest.mark.parametrize("n_variants", [100, 1000, 10000])
def test_gene_burden_analysis_scaling(
    benchmark, synthetic_gene_burden_data, gene_burden_test_config, n_variants
):
    """
    Benchmark gene burden analysis at multiple variant scales.

    Gene burden analysis aggregates variant counts per gene, performs Fisher's
    exact test, and applies multiple testing correction. This benchmark measures
    the full analysis pipeline.
    """
    # Generate synthetic gene burden data
    n_samples = 100
    n_genes = 50
    df, case_samples, control_samples = synthetic_gene_burden_data(
        n_variants, n_samples, n_genes=n_genes, seed=42
    )

    # Assign case/control counts (required columns for gene burden analysis)
    all_samples = case_samples | control_samples
    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Build config with case/control sample sets
    cfg = gene_burden_test_config.copy()
    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    # Benchmark gene burden analysis
    def run_analysis():
        return perform_gene_burden_analysis(df, cfg)

    if n_variants >= 10000:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["component"] = "gene_burden"

    # Verify result structure
    assert result is not None
    assert "GENE" in result.columns
    assert "odds_ratio" in result.columns
    assert "corrected_p_value" in result.columns


@pytest.mark.performance
@pytest.mark.parametrize("n_genes", [10, 50, 100])
def test_gene_burden_gene_scaling(
    benchmark, synthetic_gene_burden_data, gene_burden_test_config, n_genes
):
    """
    Benchmark gene burden with varying gene counts.

    Gene burden performance is affected by number of genes (multiple testing
    correction scales with gene count). This test fixes variant count and
    varies gene count.
    """
    n_variants = 1000
    n_samples = 100

    # Generate data with varying gene counts
    df, case_samples, control_samples = synthetic_gene_burden_data(
        n_variants, n_samples, n_genes=n_genes, seed=42
    )

    # Assign case/control counts
    all_samples = case_samples | control_samples
    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Build config
    cfg = gene_burden_test_config.copy()
    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    # Benchmark
    def run_analysis():
        return perform_gene_burden_analysis(df, cfg)

    if n_genes >= 100:
        result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)
    else:
        result = benchmark(run_analysis)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["component"] = "gene_burden_gene_scaling"

    # Verify
    assert result is not None
    assert len(result) <= n_genes  # Some genes may have no variants


@pytest.mark.performance
@pytest.mark.parametrize("correction_method", ["fdr", "bonferroni"])
def test_gene_burden_correction_methods(
    benchmark, synthetic_gene_burden_data, gene_burden_test_config, correction_method
):
    """
    Benchmark gene burden with different correction methods.

    Tests FDR (Benjamini-Hochberg) vs Bonferroni correction to measure
    the performance impact of different multiple testing strategies.
    """
    n_variants = 1000
    n_samples = 100
    n_genes = 50

    df, case_samples, control_samples = synthetic_gene_burden_data(
        n_variants, n_samples, n_genes=n_genes, seed=42
    )

    # Assign case/control counts
    all_samples = case_samples | control_samples
    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Build config with specific correction method
    cfg = gene_burden_test_config.copy()
    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples
    cfg["correction_method"] = correction_method

    # Benchmark
    result = benchmark(perform_gene_burden_analysis, df, cfg)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["correction_method"] = correction_method
    benchmark.extra_info["component"] = "gene_burden_correction"

    # Verify
    assert result is not None
    assert "corrected_p_value" in result.columns

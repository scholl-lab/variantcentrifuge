"""
Macro-level pipeline benchmarks.

Tests the full inheritance analysis and gene burden orchestrators at cohort scale.
These benchmarks exercise the complete analysis pipeline on realistic cohort sizes
(100-500 samples), measuring end-to-end performance of the dominant cost centers.

Macro benchmarks test Python code paths only - no external tools (bcftools, SnpSift).
The inheritance analysis orchestrator accounts for 40-60% of total pipeline time,
making it the most critical component to benchmark for optimization work.
"""

import pandas as pd
import pytest

from variantcentrifuge.gene_burden import perform_gene_burden_analysis
from variantcentrifuge.helpers import assign_case_control_counts
from variantcentrifuge.inheritance.analyzer import analyze_inheritance


def _transform_to_postprocessed_gt(df: pd.DataFrame, sample_list: list[str]) -> pd.DataFrame:
    """
    Transform GT column from comma-separated to post-replacement format.

    The inheritance analyzer expects GT column in post-genotype-replacement format:
        "SAMPLE_0001(0/1);SAMPLE_0002(0/0);SAMPLE_0003(1/1)"

    The synthetic_variants fixture generates:
        GT column with comma-separated values: "0/1,0/0,1/1,..."

    This helper transforms to the expected format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with GT column containing comma-separated genotypes
    sample_list : list[str]
        List of sample IDs in order matching GT column

    Returns
    -------
    pd.DataFrame
        DataFrame with GT column in post-replacement format
    """
    df = df.copy()

    def format_gt_row(row):
        """Convert comma-separated GT to sample(genotype) format."""
        genotypes = row["GT"].split(",")
        formatted = []
        for sample_id, gt in zip(sample_list, genotypes, strict=False):
            formatted.append(f"{sample_id}({gt})")
        return ";".join(formatted)

    df["GT"] = df.apply(format_gt_row, axis=1)
    return df


@pytest.mark.performance
def test_full_inheritance_analysis_cohort(benchmark, synthetic_variants, synthetic_pedigree):
    """
    Benchmark full inheritance analysis on medium cohort.

    Tests the complete analyze_inheritance orchestrator which performs:
    1. Per-variant pattern deduction
    2. Compound heterozygous analysis
    3. Pattern prioritization

    This is the dominant cost in the pipeline (40-60% of total runtime),
    making it the primary target for optimization.

    Cohort scale: 5000 variants x 100 samples (medium cohort)
    """
    n_variants = 5000
    n_samples = 100

    # Generate cohort pedigree (trio + 97 unrelated for realistic cohort structure)
    pedigree = synthetic_pedigree(n_samples=n_samples, seed=42)
    sample_list = list(pedigree.keys())

    # Generate synthetic data with GT column (comma-separated)
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Transform to post-replacement format (what inheritance analyzer expects)
    df = _transform_to_postprocessed_gt(df, sample_list)

    # Benchmark the full inheritance analysis
    def run_analysis():
        return analyze_inheritance(df, pedigree, sample_list, use_vectorized_comp_het=True)

    # Use pedantic mode to control iterations on expensive benchmarks
    result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)

    # Store metadata
    n_genes = df["GENE"].nunique()
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["level"] = "macro"
    benchmark.extra_info["component"] = "full_inheritance"

    # Verify result structure
    assert isinstance(result, pd.DataFrame)
    assert "Inheritance_Pattern" in result.columns
    assert len(result) == n_variants


@pytest.mark.performance
@pytest.mark.slow
def test_full_inheritance_analysis_large_cohort(benchmark, synthetic_variants, synthetic_pedigree):
    """
    Benchmark full inheritance analysis on large cohort.

    Large cohort analysis is variantcentrifuge's primary use case - rapid filtering
    of huge multi-sample VCFs for association testing and gene burden analysis.

    Cohort scale: 10000 variants x 500 samples (large cohort)

    WARNING: This benchmark may take 60+ seconds. Marked as 'slow' so it can be
    skipped with `pytest -m "performance and not slow"` during quick iterations.
    """
    n_variants = 10000
    n_samples = 500

    # Generate large cohort pedigree
    pedigree = synthetic_pedigree(n_samples=n_samples, seed=42)
    sample_list = list(pedigree.keys())

    # Generate synthetic data
    df = synthetic_variants(n_variants, n_samples, seed=42)

    # Transform to post-replacement format
    df = _transform_to_postprocessed_gt(df, sample_list)

    # Benchmark with pedantic mode (low round count due to expense)
    def run_analysis():
        return analyze_inheritance(df, pedigree, sample_list, use_vectorized_comp_het=True)

    result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)

    # Store metadata
    n_genes = df["GENE"].nunique()
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["level"] = "macro"
    benchmark.extra_info["component"] = "full_inheritance_large"

    # Verify result structure
    assert isinstance(result, pd.DataFrame)
    assert "Inheritance_Pattern" in result.columns
    assert len(result) == n_variants


@pytest.mark.performance
def test_gene_burden_full_pipeline(benchmark, synthetic_gene_burden_data, gene_burden_test_config):
    """
    Benchmark gene burden analysis full pipeline.

    Tests the complete perform_gene_burden_analysis function which performs:
    1. Per-gene variant aggregation
    2. Fisher's exact test for case/control association
    3. Multiple testing correction (FDR or Bonferroni)
    4. Odds ratio calculation

    Gene burden analysis is critical for cohort studies where the goal is
    identifying genes enriched for rare variants in cases vs controls.

    Test scale: 5000 variants x 200 samples x 100 genes
    """
    n_variants = 5000
    n_samples = 200
    n_genes = 100

    # Generate synthetic gene burden data
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

    # Use pedantic mode for controlled iteration
    result = benchmark.pedantic(run_analysis, rounds=3, iterations=1, warmup_rounds=1)

    # Store metadata
    benchmark.extra_info["n_variants"] = n_variants
    benchmark.extra_info["n_samples"] = n_samples
    benchmark.extra_info["n_genes"] = n_genes
    benchmark.extra_info["level"] = "macro"
    benchmark.extra_info["component"] = "full_gene_burden"

    # Verify result structure
    assert result is not None
    assert "GENE" in result.columns
    assert "odds_ratio" in result.columns
    assert "corrected_p_value" in result.columns

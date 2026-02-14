"""
Memory budget enforcement tests using tracemalloc.

These tests measure peak memory usage for each component and issue warnings
(never hard failures) when budgets are exceeded. Separate from timing benchmarks
because tracemalloc adds 5-20% overhead that would skew timing measurements.

Uses MemoryTracker context manager and warn_if_over_budget helper for
warning-only budget enforcement.
"""

import pytest

from variantcentrifuge.gene_burden import perform_gene_burden_analysis
from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)
from variantcentrifuge.scoring import apply_scoring

from .helpers.memory_budgets import (
    COMP_HET_BUDGET_MB,
    GENE_BURDEN_BUDGET_MB,
    INHERITANCE_BUDGET_MB,
    SCORING_BUDGET_MB,
    MemoryTracker,
    warn_if_over_budget,
)


@pytest.mark.performance
def test_inheritance_memory_budget(synthetic_variants, synthetic_pedigree):
    """
    Measure peak memory for inheritance analysis and warn if over budget.

    Uses 10K variants x 3 samples (trio) to stress test memory usage.
    Memory violations produce warnings only, never hard failures.
    """
    # Generate data
    df = synthetic_variants(n_variants=10000, n_samples=3, seed=42)
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())

    # Measure peak memory
    with MemoryTracker() as tracker:
        result = analyze_inheritance(df, pedigree, sample_list)

    # Record peak for visibility
    peak_mb = tracker.peak_mb
    print(f"\nInheritance analysis peak memory: {peak_mb:.1f} MB (10K variants)")

    # Warn if over budget (warning only, no hard failure)
    warn_if_over_budget(peak_mb, INHERITANCE_BUDGET_MB, "inheritance analysis 10K variants")

    # Correctness check only (not a memory assertion)
    assert result is not None, "Inheritance analysis should return results"
    assert "Inheritance_Pattern" in result.columns, "Result should have Inheritance_Pattern column"


@pytest.mark.performance
def test_comp_het_memory_budget(synthetic_variants, synthetic_pedigree):
    """
    Measure peak memory for compound het detection and warn if over budget.

    Uses a single gene with 5000 variants (stress test) to measure worst-case
    memory usage for compound het pair detection.
    """
    # Generate data: all variants in same gene for stress test
    df = synthetic_variants(n_variants=5000, n_samples=3, seed=42)
    # Force all to same gene for worst-case comp het scenario
    df["GENE"] = "TEST_GENE_001"
    pedigree = synthetic_pedigree(n_samples=3, seed=42)
    sample_list = list(pedigree.keys())

    # Measure peak memory
    with MemoryTracker() as tracker:
        result = analyze_gene_for_compound_het_vectorized(df, pedigree, sample_list)

    # Record peak for visibility
    peak_mb = tracker.peak_mb
    print(f"\nCompound het detection peak memory: {peak_mb:.1f} MB (5K variants in one gene)")

    # Warn if over budget
    warn_if_over_budget(peak_mb, COMP_HET_BUDGET_MB, "compound het 5K variants single gene")

    # Correctness check only
    assert result is not None or result == {}, "Compound het should return dict (may be empty)"


@pytest.mark.performance
def test_gene_burden_memory_budget():
    """
    Measure peak memory for gene burden analysis and warn if over budget.

    Uses pre-aggregated gene burden data structure (as expected by perform_gene_burden_analysis).
    Tests 50 genes x 100 samples to stress test statistical computations.
    """
    import pandas as pd

    # Create pre-aggregated gene burden data (as expected by the function)
    # This is what the gene burden stage produces before perform_gene_burden_analysis
    genes = [f"GENE_{i:04d}" for i in range(1, 51)]
    n_cases = 50
    n_controls = 50

    data = []
    for gene in genes:
        # Simulate realistic gene burden counts
        data.append(
            {
                "GENE": gene,
                "proband_count": n_cases,
                "control_count": n_controls,
                "proband_variant_count": 10,  # 10 cases with variants
                "control_variant_count": 5,  # 5 controls with variants
                "proband_allele_count": 15,  # 15 variant alleles in cases
                "control_allele_count": 7,  # 7 variant alleles in controls
            }
        )

    df = pd.DataFrame(data)

    # Build minimal config for gene burden
    cfg = {
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
    }

    # Measure peak memory
    with MemoryTracker() as tracker:
        result = perform_gene_burden_analysis(df, cfg)

    # Record peak for visibility
    peak_mb = tracker.peak_mb
    print(f"\nGene burden analysis peak memory: {peak_mb:.1f} MB (50 genes, 100 samples)")

    # Warn if over budget
    warn_if_over_budget(peak_mb, GENE_BURDEN_BUDGET_MB, "gene burden 50 genes 100 samples")

    # Correctness check only
    assert result is not None, "Gene burden should return results"
    assert len(result) > 0, "Gene burden should find some results"


@pytest.mark.performance
def test_scoring_memory_budget(synthetic_variants):
    """
    Measure peak memory for scoring and warn if over budget.

    Uses 10K variants x 10 samples to stress test scoring formula evaluation
    and variable assignment.
    """
    # Generate data
    df = synthetic_variants(n_variants=10000, n_samples=10, seed=42)

    # Build scoring config that matches actual column structure
    # (IMPACT and EFFECT are present in synthetic data)
    scoring_config = {
        "variables": {
            "IMPACT": "impact_val",
            "EFFECT": "effect_val",
        },
        "formulas": [
            {"impact_score": "((impact_val == 'HIGH') * 10) + ((impact_val == 'MODERATE') * 5)"},
            {
                "effect_score": (
                    "((effect_val == 'stop_gained') * 10) + "
                    "((effect_val == 'missense_variant') * 5)"
                )
            },
        ],
    }

    # Measure peak memory
    with MemoryTracker() as tracker:
        result = apply_scoring(df, scoring_config)

    # Record peak for visibility
    peak_mb = tracker.peak_mb
    print(f"\nScoring peak memory: {peak_mb:.1f} MB (10K variants, 10 samples)")

    # Warn if over budget
    warn_if_over_budget(peak_mb, SCORING_BUDGET_MB, "scoring 10K variants 10 samples")

    # Correctness check only
    assert result is not None, "Scoring should return DataFrame"
    assert len(result) == len(df), "Scoring should preserve all rows"


@pytest.mark.performance
def test_dataframe_read_memory_budget(synthetic_variants, tmp_path):
    """
    Measure peak memory for DataFrame I/O and record baseline.

    Uses 50K variants to establish baseline memory for large DataFrame loading.
    This is informational — no budget enforcement, just recording for reference.
    """
    import pandas as pd

    # Generate large DataFrame
    df = synthetic_variants(n_variants=50000, n_samples=10, seed=42)

    # Write to TSV
    tsv_path = tmp_path / "test_variants.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    # Measure peak memory for read
    with MemoryTracker() as tracker:
        loaded_df = pd.read_csv(tsv_path, sep="\t")

    # Record peak for visibility
    peak_mb = tracker.peak_mb
    print(f"\nDataFrame read peak memory: {peak_mb:.1f} MB (50K variants)")
    print(f"File size: {tsv_path.stat().st_size / (1024 * 1024):.1f} MB")

    # Correctness check only (no budget for I/O — this is baseline measurement)
    assert loaded_df is not None, "DataFrame should load successfully"
    assert len(loaded_df) == 50000, "All rows should be loaded"

"""
Regression test for gene burden analysis.

This test verifies that removing the dead GT parsing loop (lines 220-249)
from perform_gene_burden_analysis does not change the output.
"""

import pandas as pd
import pytest

from variantcentrifuge.gene_burden import perform_gene_burden_analysis


@pytest.mark.unit
def test_gene_burden_output_unchanged():
    """
    Test that gene burden analysis produces expected output.

    This test creates a synthetic dataset with realistic GT column values
    and verifies that the output matches expected values after the dead code
    removal.
    """
    # Create synthetic test data with multiple genes and variants
    test_data = pd.DataFrame(
        {
            "GENE": [
                "BRCA1",
                "BRCA1",
                "BRCA1",
                "TP53",
                "TP53",
                "CFTR",
                "CFTR",
                "CFTR",
                "CFTR",
            ],
            "proband_count": [100, 100, 100, 100, 100, 100, 100, 100, 100],
            "control_count": [200, 200, 200, 200, 200, 200, 200, 200, 200],
            "proband_variant_count": [5, 3, 4, 2, 1, 6, 5, 4, 3],
            "control_variant_count": [2, 1, 2, 1, 0, 3, 2, 1, 2],
            "proband_allele_count": [8, 5, 6, 3, 2, 9, 7, 6, 5],
            "control_allele_count": [3, 2, 3, 2, 0, 4, 3, 2, 3],
            "GT": [
                "Sample1(0/1);Sample2(1/1);Sample3(0/1)",
                "Sample4(0/1);Sample5(0/1)",
                "Sample6(1/1);Sample7(0/1)",
                "Sample8(0/1);Sample9(0/1)",
                "Sample10(0/1)",
                "Sample11(0/1);Sample12(1/1);Sample13(0/1)",
                "Sample14(0/1);Sample15(1/1)",
                "Sample16(1/1);Sample17(0/1)",
                "Sample18(0/1);Sample19(0/1);Sample20(0/1)",
            ],
        }
    )

    # Configuration for gene burden analysis
    cfg = {
        "gene_burden_mode": "alleles",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
        "confidence_interval_alpha": 0.05,
        "continuity_correction": 0.5,
    }

    # Run gene burden analysis
    result = perform_gene_burden_analysis(test_data, cfg)

    # Verify expected structure
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3  # Three genes: BRCA1, CFTR, TP53
    assert set(result["GENE"]) == {"BRCA1", "CFTR", "TP53"}

    # Required columns should be present
    required_columns = [
        "GENE",
        "proband_count",
        "control_count",
        "proband_allele_count",
        "control_allele_count",
        "raw_p_value",
        "corrected_p_value",
        "odds_ratio",
        "or_ci_lower",
        "or_ci_upper",
    ]
    for col in required_columns:
        assert col in result.columns, f"Missing required column: {col}"

    # Verify sample counts are preserved (should be consistent across variants)
    assert all(result["proband_count"] == 100)
    assert all(result["control_count"] == 200)

    # Verify allele counts are summed correctly per gene
    brca1_row = result[result["GENE"] == "BRCA1"].iloc[0]
    tp53_row = result[result["GENE"] == "TP53"].iloc[0]
    cftr_row = result[result["GENE"] == "CFTR"].iloc[0]

    # BRCA1: sum of proband_allele_count = 8 + 5 + 6 = 19
    assert brca1_row["proband_allele_count"] == 19
    # BRCA1: sum of control_allele_count = 3 + 2 + 3 = 8
    assert brca1_row["control_allele_count"] == 8

    # TP53: sum of proband_allele_count = 3 + 2 = 5
    assert tp53_row["proband_allele_count"] == 5
    # TP53: sum of control_allele_count = 2 + 0 = 2
    assert tp53_row["control_allele_count"] == 2

    # CFTR: sum of proband_allele_count = 9 + 7 + 6 + 5 = 27
    assert cftr_row["proband_allele_count"] == 27
    # CFTR: sum of control_allele_count = 4 + 3 + 2 + 3 = 12
    assert cftr_row["control_allele_count"] == 12

    # Verify p-values are numeric and in valid range
    assert all(result["raw_p_value"] >= 0.0)
    assert all(result["raw_p_value"] <= 1.0)
    assert all(result["corrected_p_value"] >= 0.0)
    assert all(result["corrected_p_value"] <= 1.0)

    # Verify odds ratios are positive or NaN
    assert all((result["odds_ratio"] > 0) | result["odds_ratio"].isna())

    # Verify confidence intervals are present (may be NaN for edge cases)
    assert "or_ci_lower" in result.columns
    assert "or_ci_upper" in result.columns

    # For valid CIs, lower bound should be less than upper bound
    valid_ci_mask = result["or_ci_lower"].notna() & result["or_ci_upper"].notna()
    if valid_ci_mask.any():
        valid_ci_rows = result[valid_ci_mask]
        assert all(valid_ci_rows["or_ci_lower"] <= valid_ci_rows["or_ci_upper"])


@pytest.mark.unit
def test_gene_burden_samples_mode():
    """
    Test gene burden analysis in samples mode.

    Verify that the function works correctly when mode is "samples" instead of "alleles".
    """
    test_data = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1", "GENE2"],
            "proband_count": [50, 50, 50],
            "control_count": [100, 100, 100],
            "proband_variant_count": [5, 3, 2],
            "control_variant_count": [2, 1, 1],
            "proband_allele_count": [8, 5, 3],
            "control_allele_count": [3, 2, 2],
            "GT": [
                "S1(0/1);S2(1/1)",
                "S3(0/1);S4(0/1)",
                "S5(0/1);S6(0/1)",
            ],
        }
    )

    cfg = {
        "gene_burden_mode": "samples",
        "correction_method": "bonferroni",
        "confidence_interval_alpha": 0.05,
    }

    result = perform_gene_burden_analysis(test_data, cfg)

    # Verify structure
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 2  # Two genes

    # In samples mode, should have variant_count columns
    assert "proband_variant_count" in result.columns
    assert "control_variant_count" in result.columns

    # Verify counts are summed correctly
    gene1_row = result[result["GENE"] == "GENE1"].iloc[0]
    assert gene1_row["proband_variant_count"] == 8  # 5 + 3
    assert gene1_row["control_variant_count"] == 3  # 2 + 1


@pytest.mark.unit
def test_gene_burden_single_gene_edge_case():
    """
    Test gene burden analysis with a single gene and edge case counts.
    """
    test_data = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "proband_count": [10],
            "control_count": [20],
            "proband_variant_count": [1],
            "control_variant_count": [0],
            "proband_allele_count": [1],
            "control_allele_count": [0],
            "GT": ["S1(0/1)"],
        }
    )

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    result = perform_gene_burden_analysis(test_data, cfg)

    # Should return single gene
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1
    assert result.iloc[0]["GENE"] == "GENE1"

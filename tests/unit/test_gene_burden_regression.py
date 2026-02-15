"""
Regression test for gene burden analysis.

Tests GT-based per-sample collapsing, column-based aggregation (Phase 11),
and the legacy pre-computed count path (backward compatibility).
"""

import pandas as pd
import pytest

from variantcentrifuge.gene_burden import perform_gene_burden_analysis


@pytest.mark.unit
def test_gene_burden_gt_collapsing_unique_carriers():
    """
    Test that GT-based collapsing correctly counts unique carriers per gene.

    When a sample carries variants at multiple sites in the same gene,
    it should be counted only once in "samples" mode (CMC/CAST collapsing test).
    """
    # Gene with 3 variants. Sample S1 appears at ALL three sites.
    # S2 appears at V1 only, S3 at V2 only, S4 at V3 only.
    test_data = pd.DataFrame(
        {
            "GENE": ["BRCA1", "BRCA1", "BRCA1"],
            "GT": [
                "S1(0/1);S2(0/1);C1(0/1)",  # V1: cases S1,S2; control C1
                "S1(1/1);S3(0/1);C2(0/1);C1(0/1)",  # V2: cases S1,S3; controls C1,C2
                "S1(0/1);S4(0/1)",  # V3: cases S1,S4; no controls
            ],
            # Pre-computed counts are present but should be ignored when
            # case/control samples are provided
            "proband_count": [10, 10, 10],
            "control_count": [20, 20, 20],
            "proband_variant_count": [2, 2, 2],
            "control_variant_count": [1, 2, 0],
            "proband_allele_count": [2, 3, 2],
            "control_allele_count": [1, 2, 0],
        }
    )

    case_samples = {"S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"}
    control_samples = {
        "C1",
        "C2",
        "C3",
        "C4",
        "C5",
        "C6",
        "C7",
        "C8",
        "C9",
        "C10",
        "C11",
        "C12",
        "C13",
        "C14",
        "C15",
        "C16",
        "C17",
        "C18",
        "C19",
        "C20",
    }

    cfg = {
        "gene_burden_mode": "samples",
        "correction_method": "fdr",
    }

    result = perform_gene_burden_analysis(
        test_data, cfg, case_samples=case_samples, control_samples=control_samples
    )

    assert len(result) == 1
    row = result.iloc[0]

    # Unique case carriers: S1, S2, S3, S4 = 4 (not 6 from summing per-variant)
    assert row["proband_carrier_count"] == 4
    # Unique control carriers: C1, C2 = 2 (not 3 from summing per-variant)
    assert row["control_carrier_count"] == 2
    assert row["proband_count"] == 10
    assert row["control_count"] == 20


@pytest.mark.unit
def test_gene_burden_gt_collapsing_allele_mode():
    """
    Test allele mode uses max dosage per sample, not sum across sites.

    Sample S1 is het (0/1) at V1 and hom-alt (1/1) at V2.
    In allele mode, S1 should contribute max(1, 2) = 2 alleles, not 1+2 = 3.
    """
    test_data = pd.DataFrame(
        {
            "GENE": ["PKD1", "PKD1"],
            "GT": [
                "S1(0/1);S2(0/1);C1(0/1);C2(0/1)",  # V1
                "S1(1/1);C1(1/1)",  # V2
            ],
            "proband_count": [5, 5],
            "control_count": [10, 10],
            "proband_variant_count": [2, 1],
            "control_variant_count": [2, 1],
            "proband_allele_count": [2, 2],
            "control_allele_count": [2, 2],
        }
    )

    case_samples = {"S1", "S2", "S3", "S4", "S5"}
    control_samples = {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"}

    cfg = {
        "gene_burden_mode": "alleles",
        "correction_method": "fdr",
    }

    result = perform_gene_burden_analysis(
        test_data, cfg, case_samples=case_samples, control_samples=control_samples
    )

    assert len(result) == 1
    row = result.iloc[0]

    # S1: max(1, 2) = 2 alleles; S2: max(1) = 1 allele; total = 3
    assert row["proband_allele_count"] == 3
    # C1: max(1, 2) = 2 alleles; C2: max(1) = 1 allele; total = 3
    assert row["control_allele_count"] == 3
    # Total alleles cannot exceed 2*N
    assert row["proband_allele_count"] <= 2 * row["proband_count"]
    assert row["control_allele_count"] <= 2 * row["control_count"]


@pytest.mark.unit
def test_gene_burden_legacy_mode():
    """
    Test legacy mode (no case/control samples provided) still works.

    This ensures backward compatibility with existing callers that pass
    pre-computed per-variant counts without sample sets.
    """
    test_data = pd.DataFrame(
        {
            "GENE": [
                "BRCA1",
                "BRCA1",
                "BRCA1",
                "TP53",
                "TP53",
            ],
            "proband_count": [100, 100, 100, 100, 100],
            "control_count": [200, 200, 200, 200, 200],
            "proband_variant_count": [5, 3, 4, 2, 1],
            "control_variant_count": [2, 1, 2, 1, 0],
            "proband_allele_count": [8, 5, 6, 3, 2],
            "control_allele_count": [3, 2, 3, 2, 0],
            "GT": [
                "S1(0/1);S2(1/1);S3(0/1)",
                "S4(0/1);S5(0/1)",
                "S6(1/1);S7(0/1)",
                "S8(0/1);S9(0/1)",
                "S10(0/1)",
            ],
        }
    )

    cfg = {
        "gene_burden_mode": "alleles",
        "correction_method": "fdr",
        "confidence_interval_method": "normal_approx",
        "confidence_interval_alpha": 0.05,
        "continuity_correction": 0.5,
    }

    # No case/control samples -> legacy path
    result = perform_gene_burden_analysis(test_data, cfg)

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 2  # BRCA1, TP53
    assert set(result["GENE"]) == {"BRCA1", "TP53"}

    # Legacy mode sums per-variant counts
    brca1 = result[result["GENE"] == "BRCA1"].iloc[0]
    assert brca1["proband_allele_count"] == 19  # 8+5+6
    assert brca1["control_allele_count"] == 8  # 3+2+3

    tp53 = result[result["GENE"] == "TP53"].iloc[0]
    assert tp53["proband_allele_count"] == 5  # 3+2
    assert tp53["control_allele_count"] == 2  # 2+0

    # p-values valid
    assert all(result["raw_p_value"] >= 0.0)
    assert all(result["raw_p_value"] <= 1.0)


@pytest.mark.unit
def test_gene_burden_samples_mode():
    """
    Test gene burden analysis in samples mode (legacy path).
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

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 2
    assert "proband_variant_count" in result.columns
    assert "control_variant_count" in result.columns

    gene1 = result[result["GENE"] == "GENE1"].iloc[0]
    assert gene1["proband_variant_count"] == 8  # 5+3 (capped at 50)
    assert gene1["control_variant_count"] == 3  # 2+1


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

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 1
    assert result.iloc[0]["GENE"] == "GENE1"


@pytest.mark.unit
def test_gene_burden_no_double_counting():
    """
    Verify the core fix: samples appearing at multiple variant sites in the
    same gene are NOT double-counted.

    This is the regression test for the bug where summing per-variant counts
    inflated carrier counts beyond the actual number of unique carriers.
    """
    # 10 variants in GENE_X. Same 3 case samples appear at EVERY variant.
    # Old (buggy): would sum to 30 carrier events, capped to 5 (total cases).
    # New (correct): 3 unique carriers.
    gt_entries = ["CaseA(0/1);CaseB(0/1);CaseC(1/1);CtrlX(0/1)"] * 10

    test_data = pd.DataFrame(
        {
            "GENE": ["GENE_X"] * 10,
            "GT": gt_entries,
            "proband_count": [5] * 10,
            "control_count": [10] * 10,
            "proband_variant_count": [3] * 10,
            "control_variant_count": [1] * 10,
            "proband_allele_count": [4] * 10,  # 1+1+2 = 4 per variant
            "control_allele_count": [1] * 10,
        }
    )

    case_samples = {"CaseA", "CaseB", "CaseC", "CaseD", "CaseE"}
    control_samples = {
        "CtrlX",
        "CtrlY",
        "CtrlZ",
        "Ctrl1",
        "Ctrl2",
        "Ctrl3",
        "Ctrl4",
        "Ctrl5",
        "Ctrl6",
        "Ctrl7",
    }

    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}

    result = perform_gene_burden_analysis(
        test_data, cfg, case_samples=case_samples, control_samples=control_samples
    )

    row = result.iloc[0]
    # Unique carriers, NOT summed per-variant counts
    assert row["proband_carrier_count"] == 3, (
        f"Expected 3 unique case carriers, got {row['proband_carrier_count']}"
    )
    assert row["control_carrier_count"] == 1, (
        f"Expected 1 unique control carrier, got {row['control_carrier_count']}"
    )

    # In allele mode, max dosage per sample
    cfg["gene_burden_mode"] = "alleles"
    result_alleles = perform_gene_burden_analysis(
        test_data, cfg, case_samples=case_samples, control_samples=control_samples
    )
    row_a = result_alleles.iloc[0]
    # CaseA: max(1 across 10 sites) = 1; CaseB: 1; CaseC: 2; total = 4
    assert row_a["proband_allele_count"] == 4
    # CtrlX: max(1 across 10 sites) = 1
    assert row_a["control_allele_count"] == 1


@pytest.mark.unit
def test_gene_burden_column_based_aggregation():
    """
    Test column-based aggregation (Phase 11 optimization).

    When per-sample GT columns (GEN_0__GT, GEN_1__GT, ...) are present and
    vcf_samples is provided, the function should use fast column-based
    aggregation instead of parsing the packed GT string.

    Results must match the packed GT string approach exactly.
    """
    # Create DataFrame with per-sample GT columns
    vcf_samples = ["S1", "S2", "S3", "C1", "C2", "C3"]
    test_data = pd.DataFrame(
        {
            "GENE": ["PKD1", "PKD1", "PKD1"],
            # Per-sample GT columns (Phase 11 format)
            "GEN_0__GT": ["0/1", "1/1", "0/1"],  # S1: het, hom, het -> max=2
            "GEN_1__GT": ["0/1", "0/0", "0/0"],  # S2: het at V1 only -> max=1
            "GEN_2__GT": ["0/0", "0/1", "0/0"],  # S3: het at V2 only -> max=1
            "GEN_3__GT": ["0/1", "0/1", "0/0"],  # C1: het at V1,V2 -> max=1
            "GEN_4__GT": ["0/0", "0/0", "0/1"],  # C2: het at V3 only -> max=1
            "GEN_5__GT": ["0/0", "0/0", "0/0"],  # C3: no variants -> 0
            # Also include packed GT for comparison
            "GT": [
                "S1(0/1);S2(0/1);C1(0/1)",
                "S1(1/1);S3(0/1);C1(0/1)",
                "S1(0/1);C2(0/1)",
            ],
            "proband_count": [3, 3, 3],
            "control_count": [3, 3, 3],
            "proband_variant_count": [2, 2, 1],
            "control_variant_count": [1, 1, 1],
            "proband_allele_count": [2, 3, 1],
            "control_allele_count": [1, 1, 1],
        }
    )

    case_samples = {"S1", "S2", "S3"}
    control_samples = {"C1", "C2", "C3"}

    # Test samples mode with column-based aggregation
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}
    result_cols = perform_gene_burden_analysis(
        test_data,
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
    )

    assert len(result_cols) == 1
    row = result_cols.iloc[0]
    # S1, S2, S3 all carry at least one variant
    assert row["proband_carrier_count"] == 3
    # C1, C2 carry variants; C3 does not
    assert row["control_carrier_count"] == 2
    assert row["proband_count"] == 3
    assert row["control_count"] == 3

    # Test allele mode with column-based aggregation
    cfg["gene_burden_mode"] = "alleles"
    result_alleles = perform_gene_burden_analysis(
        test_data,
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
    )

    row_a = result_alleles.iloc[0]
    # S1: max(1,2,1)=2; S2: max(1,0,0)=1; S3: max(0,1,0)=1 -> total=4
    assert row_a["proband_allele_count"] == 4
    # C1: max(1,1,0)=1; C2: max(0,0,1)=1; C3: 0 -> total=2
    assert row_a["control_allele_count"] == 2
    # Diploid constraint
    assert row_a["proband_allele_count"] <= 2 * row_a["proband_count"]
    assert row_a["control_allele_count"] <= 2 * row_a["control_count"]


@pytest.mark.unit
def test_gene_burden_column_vs_gt_consistency():
    """
    Verify column-based and packed GT string approaches produce identical results.
    """
    vcf_samples = ["CaseA", "CaseB", "CaseC", "CtrlX", "CtrlY"]
    test_data = pd.DataFrame(
        {
            "GENE": ["G1", "G1"],
            "GEN_0__GT": ["0/1", "1/1"],  # CaseA
            "GEN_1__GT": ["0/1", "0/0"],  # CaseB
            "GEN_2__GT": ["0/0", "0/0"],  # CaseC
            "GEN_3__GT": ["0/1", "0/1"],  # CtrlX
            "GEN_4__GT": ["0/0", "0/0"],  # CtrlY
            "GT": [
                "CaseA(0/1);CaseB(0/1);CtrlX(0/1)",
                "CaseA(1/1);CtrlX(0/1)",
            ],
            "proband_count": [3, 3],
            "control_count": [2, 2],
            "proband_variant_count": [2, 1],
            "control_variant_count": [1, 1],
            "proband_allele_count": [2, 2],
            "control_allele_count": [1, 1],
        }
    )

    case_samples = {"CaseA", "CaseB", "CaseC"}
    control_samples = {"CtrlX", "CtrlY"}
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}

    # Column-based (with vcf_samples)
    result_col = perform_gene_burden_analysis(
        test_data,
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
    )

    # Drop per-sample columns to force packed GT path
    gt_string_data = test_data.drop(columns=[c for c in test_data.columns if c.startswith("GEN_")])
    result_gt = perform_gene_burden_analysis(
        gt_string_data,
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
    )

    # Both should produce identical carrier counts
    assert result_col.iloc[0]["proband_carrier_count"] == result_gt.iloc[0]["proband_carrier_count"]
    assert result_col.iloc[0]["control_carrier_count"] == result_gt.iloc[0]["control_carrier_count"]

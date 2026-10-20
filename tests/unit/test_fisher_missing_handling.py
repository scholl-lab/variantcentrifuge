"""
Unit tests for Fisher exact test missing genotype handling (#138).

Verifies that missing genotypes (./., ., .|., None, pd.NA, NA) are:
1. Preserved during GT reconstruction rather than omitted/skipped.
2. Excluded from reference allele/carrier counts in Fisher exact contingency tables
   under a complete-case policy across both gene_burden.py and association/tests/fisher.py.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.base import AssociationConfig
from variantcentrifuge.association.tests.fisher import FisherExactTest
from variantcentrifuge.gene_burden import (
    _aggregate_gene_burden_from_columns,
    _aggregate_gene_burden_from_gt,
    perform_gene_burden_analysis,
)
from variantcentrifuge.gt_reconstruction import reconstruct_gt_column


@pytest.mark.unit
def test_reconstruct_gt_preserves_missing_calls():
    """Verify reconstruct_gt_column preserves missing genotypes and elides reference."""
    df = pd.DataFrame(
        {
            "CHROM": ["chr1", "chr1", "chr1"],
            "POS": [100, 200, 300],
            "GEN_0__GT": ["0/1", "1/1", ""],
            "GEN_1__GT": [None, pd.NA, "."],
            "GEN_2__GT": ["1/0", "0/0", "0/1"],
            "GEN_3__GT": ["./.", ".|.", "NA"],
        }
    )
    vcf_samples = ["S1", "S2", "S3", "S4"]

    res = reconstruct_gt_column(df, vcf_samples)

    assert "GT" in res.columns
    # Row 0: S1(0/1), S2(./.), S3(1/0), S4(./.)
    assert res["GT"].iloc[0] == "S1(0/1);S2(./.);S3(1/0);S4(./.)"
    # Row 1: S1(1/1), S2(./.), S3(0/0 is elided), S4(.|.)
    assert res["GT"].iloc[1] == "S1(1/1);S2(./.);S4(.|.)"
    # Row 2: S1(./.), S2(.), S3(0/1), S4(./.)
    assert res["GT"].iloc[2] == "S1(./.);S2(.);S3(0/1);S4(./.)"


@pytest.mark.unit
def test_fisher_missing_handling_single_variant():
    """Cases [0/1, ./.], controls [0/1, 0/0] must produce table [[1, 1], [1, 3]] and OR=3."""
    # Column-based and packed-GT test
    df = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "CHROM": ["chr1"],
            "POS": [1000],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["./."],
            "GEN_2__GT": ["0/1"],
            "GEN_3__GT": ["0/0"],
            "GT": ["Case1(0/1);Case2(./.);Ctrl1(0/1)"],
        }
    )
    case_samples = {"Case1", "Case2"}
    control_samples = {"Ctrl1", "Ctrl2"}
    vcf_samples = ["Case1", "Case2", "Ctrl1", "Ctrl2"]

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    # 1. perform_gene_burden_analysis (column-based)
    res_col = perform_gene_burden_analysis(
        df.copy(),
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
    )
    assert len(res_col) == 1
    row_col = res_col.iloc[0]
    assert row_col["proband_allele_count"] == 1
    assert row_col["control_allele_count"] == 1
    assert row_col["proband_count"] == 1  # 1 called case
    assert row_col["control_count"] == 2  # 2 called controls
    assert np.isclose(row_col["odds_ratio"], 3.0)

    # 2. perform_gene_burden_analysis (packed GT string fallback)
    df_no_gen = df[["GENE", "CHROM", "POS", "GT"]].copy()
    res_gt = perform_gene_burden_analysis(
        df_no_gen, cfg, case_samples=case_samples, control_samples=control_samples
    )
    assert len(res_gt) == 1
    row_gt = res_gt.iloc[0]
    assert row_gt["proband_allele_count"] == 1
    assert row_gt["control_allele_count"] == 1
    assert row_gt["proband_count"] == 1
    assert row_gt["control_count"] == 2
    assert np.isclose(row_gt["odds_ratio"], 3.0)

    # 3. FisherExactTest from association test engine using actual aggregator outputs
    gt_columns = ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"]
    agg_cols = _aggregate_gene_burden_from_columns(
        df, case_samples, control_samples, vcf_samples, gt_columns
    )[0]
    agg_gt = _aggregate_gene_burden_from_gt(df_no_gen, case_samples, control_samples)[0]

    assert agg_cols["proband_count"] == 2
    assert agg_cols["proband_called_count"] == 1
    assert agg_cols["proband_called_chromosomes"] == 2
    assert agg_cols["control_called_chromosomes"] == 4

    assert agg_gt["proband_count"] == 2
    assert agg_gt["proband_called_count"] == 1
    assert agg_gt["proband_called_chromosomes"] == 2
    assert agg_gt["control_called_chromosomes"] == 4

    assoc_cfg = AssociationConfig(gene_burden_mode="alleles")
    fisher_test = FisherExactTest()

    res_col = fisher_test.run("GENE1", agg_cols, assoc_cfg)
    assert res_col.effect_size is not None
    assert np.isclose(res_col.effect_size, 3.0)
    assert res_col.extra["table"] == [[1, 1], [1, 3]]

    res_gt_fisher = fisher_test.run("GENE1", agg_gt, assoc_cfg)
    assert res_gt_fisher.effect_size is not None
    assert np.isclose(res_gt_fisher.effect_size, 3.0)
    assert res_gt_fisher.extra["table"] == [[1, 1], [1, 3]]

    # 4. Haploid denominator assertion with FisherExactTest on aggregator output
    df_hap = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "CHROM": ["chrX"],
            "POS": [1000],
            "GEN_0__GT": ["1"],
            "GEN_1__GT": ["./."],
            "GEN_2__GT": ["0/1"],
            "GEN_3__GT": ["0/0"],
        }
    )
    agg_hap = _aggregate_gene_burden_from_columns(
        df_hap, case_samples, control_samples, vcf_samples, gt_columns
    )[0]
    assert agg_hap["proband_called_chromosomes"] == 1
    assert agg_hap["control_called_chromosomes"] == 4
    res_hap = fisher_test.run("GENE1", agg_hap, assoc_cfg)
    assert res_hap.extra["table"] == [[1, 1], [0, 3]]


@pytest.mark.unit
def test_fisher_missing_handling_multi_variant_complete_case():
    """Multi-variant fixture: sample with ./. is excluded under complete-case."""
    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "CHROM": ["chr1", "chr1"],
            "POS": [1000, 2000],
            "GEN_0__GT": ["0/1", "0/0"],  # Case1: complete carrier (alt at V1, ref at V2)
            "GEN_1__GT": ["0/0", "./."],  # Case2: incomplete (missing at V2) -> excluded
            "GEN_2__GT": ["0/1", "0/0"],  # Ctrl1: complete carrier
            "GEN_3__GT": ["0/0", "0/0"],  # Ctrl2: complete ref
            "GT": [
                "Case1(0/1);Ctrl1(0/1)",
                "Case2(./.)",
            ],
        }
    )
    case_samples = {"Case1", "Case2"}
    control_samples = {"Ctrl1", "Ctrl2"}
    vcf_samples = ["Case1", "Case2", "Ctrl1", "Ctrl2"]

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    res = perform_gene_burden_analysis(
        df, cfg, case_samples=case_samples, control_samples=control_samples, vcf_samples=vcf_samples
    )
    assert len(res) == 1
    row = res.iloc[0]
    assert row["proband_allele_count"] == 1
    assert row["control_allele_count"] == 1
    assert row["proband_count"] == 1
    assert row["control_count"] == 2
    assert np.isclose(row["odds_ratio"], 3.0)


@pytest.mark.unit
def test_fisher_all_missing_handling():
    """All cases missing should skip in gene_burden and return p_value=None in FisherExactTest."""
    df = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "CHROM": ["chr1"],
            "POS": [1000],
            "GEN_0__GT": ["./."],
            "GEN_1__GT": ["./."],
            "GEN_2__GT": ["0/1"],
            "GEN_3__GT": ["0/0"],
            "GT": ["Case1(./.);Case2(./.);Ctrl1(0/1)"],
        }
    )
    case_samples = {"Case1", "Case2"}
    control_samples = {"Ctrl1", "Ctrl2"}
    vcf_samples = ["Case1", "Case2", "Ctrl1", "Ctrl2"]

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    res = perform_gene_burden_analysis(
        df, cfg, case_samples=case_samples, control_samples=control_samples, vcf_samples=vcf_samples
    )
    # When all cases are missing (p_count == 0), the gene is not testable and skipped
    assert res.empty

    # FisherExactTest returns p_value=None when proband_count == 0 or control_count == 0
    fisher_test = FisherExactTest()
    assoc_cfg = AssociationConfig(gene_burden_mode="alleles")

    res_cases_zero = fisher_test.run(
        "GENE1",
        {
            "GENE": "GENE1",
            "proband_count": 0,
            "control_count": 2,
            "proband_allele_count": 0,
            "control_allele_count": 1,
            "n_qualifying_variants": 1,
        },
        assoc_cfg,
    )
    assert res_cases_zero.p_value is None

    res_ctrls_zero = fisher_test.run(
        "GENE1",
        {
            "GENE": "GENE1",
            "proband_count": 2,
            "control_count": 0,
            "proband_allele_count": 1,
            "control_allele_count": 0,
            "n_qualifying_variants": 1,
        },
        assoc_cfg,
    )
    assert res_ctrls_zero.p_value is None


@pytest.mark.unit
def test_reconstruct_gt_float_nan_and_phased_ref():
    """Verify float('nan') and 0|0 in reconstruction and collapsing."""
    from variantcentrifuge.gt_reconstruction import (
        _collapse_per_sample_gt_columns_for_output,
        reconstruct_gt_column,
    )

    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [100],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": [float("nan")],
            "GEN_2__GT": ["0|0"],
            "GEN_3__GT": ["1/1"],
        }
    )
    vcf_samples = ["S1", "S2", "S3", "S4"]

    res = reconstruct_gt_column(df.copy(), vcf_samples)
    # S1(0/1) included, S2(nan -> ./.) included, S3(0|0) elided, S4(1/1) included
    assert res["GT"].iloc[0] == "S1(0/1);S2(./.);S4(1/1)"

    # Test _collapse_per_sample_gt_columns_for_output with pseudonym mapping
    df_with_gt = df.copy()
    df_with_gt["GT"] = "P1(0/1);P2(./.);P4(1/1)"
    collapsed = _collapse_per_sample_gt_columns_for_output(
        df_with_gt,
        vcf_samples=vcf_samples,
    )
    assert "GEN_0__GT" not in collapsed.columns
    assert collapsed["GT"].iloc[0] == "P1(0/1);P2(./.);P4(1/1)"


@pytest.mark.unit
def test_phenotype_map_does_not_inflate_carriers():
    """Verify missing and ref calls do not assign phenotypes to samples."""
    from variantcentrifuge.phenotype_helpers import build_sample_phenotype_map

    df = pd.DataFrame(
        {
            "GT": ["S1(0/1);S2(./.);S3(0/0);S4(0|0);S5(1/1)"],
            "phenotypes": ["phenoA;phenoB;phenoC;phenoD;phenoE"],
        }
    )

    sample_phenos = build_sample_phenotype_map(df)
    # Only S1 (0/1) and S5 (1/1) are variant calls; S2(./.), S3(0/0), S4(0|0) must be skipped.
    # Since only 2 samples are variants but 5 pheno groups are present, counts mismatch (2 != 5),
    # so phenotype assignment must not wrongly align missing/ref samples.
    assert "S2" not in sample_phenos
    assert "S3" not in sample_phenos
    assert "S4" not in sample_phenos


@pytest.mark.unit
def test_ploidy_aware_chromosome_counts_in_allele_mode():
    """Verify haploid calls (1 chrom) vs diploid calls (2 chroms) in allele mode."""
    # Gene with 1 haploid case (genotype '1') and 1 diploid control (genotype '0/0')
    df = pd.DataFrame(
        {
            "GENE": ["CHRX_GENE"],
            "CHROM": ["chrX"],
            "POS": [1000],
            "GEN_0__GT": ["1"],
            "GEN_1__GT": ["0/0"],
        }
    )
    case_samples = {"Case1"}
    control_samples = {"Ctrl1"}
    vcf_samples = ["Case1", "Ctrl1"]

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    res = perform_gene_burden_analysis(
        df, cfg, case_samples=case_samples, control_samples=control_samples, vcf_samples=vcf_samples
    )
    assert len(res) == 1
    row = res.iloc[0]
    assert row["proband_allele_count"] == 1
    assert row["control_allele_count"] == 0
    # Case has 1 chromosome called, 1 alt allele -> 0 ref alleles
    # Control has 2 chromosomes called, 0 alt allele -> 2 ref alleles
    # Table: [[1, 0], [0, 2]] -> OR is infinity or nan depending on continuity
    assert not np.isnan(row["raw_p_value"])


@pytest.mark.unit
def test_legacy_aggregation_count_only_multi_site_missing_rejected(caplog):
    """Verify count-only multi-site missingness without sample IDs is skipped (p_called=0)."""
    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "proband_count": [2, 2],
            "control_count": [2, 2],
            "proband_called_count": [1, 1],
            "control_called_count": [2, 2],
            "proband_variant_count": [1, 1],
            "control_variant_count": [0, 0],
            "proband_allele_count": [1, 1],
            "control_allele_count": [0, 0],
        }
    )
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}
    with caplog.at_level("WARNING", logger="variantcentrifuge"):
        res = perform_gene_burden_analysis(df, cfg)
    # Untestable multi-site count-only missingness is rejected -> gene skipped, loudly
    assert res.empty
    assert "Gene GENE1 skipped: legacy counts cannot tell which samples are missing" in caplog.text


@pytest.mark.unit
def test_legacy_aggregation_count_only_single_site_missing_allowed():
    """Verify count-only single-site gene missingness uses minimum called count."""
    df = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "proband_count": [2],
            "control_count": [2],
            "proband_called_count": [1],
            "control_called_count": [2],
            "proband_variant_count": [1],
            "control_variant_count": [0],
            "proband_allele_count": [1],
            "control_allele_count": [0],
        }
    )
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}
    res = perform_gene_burden_analysis(df, cfg)
    assert len(res) == 1
    row = res.iloc[0]
    assert row["proband_count"] == 1
    assert row["control_count"] == 2


@pytest.mark.unit
def test_legacy_aggregation_with_missing_samples_tracking():
    """Verify legacy aggregation via assign_case_control_counts tracks sample IDs across sites.

    - Same sample missing at 2 sites: union = {P1}, retains 2 - 1 = 1 case.
    - 2 different samples missing at 2 sites: union = {P1, P2}, retains 2 - 2 = 0 cases (skipped).
    """
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts

    case_samples = {"P1", "P2"}
    control_samples = {"C1", "C2"}
    all_samples = {"P1", "P2", "C1", "C2"}
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}

    # Case A: Same sample (P1) missing at both sites -> 1 case retained
    df_same = pd.DataFrame(
        {
            "GENE": ["GENE_SAME", "GENE_SAME"],
            "GT": [
                "P1(./.);P2(0/1);C1(0/0);C2(0/0)",
                "P1(./.);P2(0/0);C1(0/0);C2(0/0)",
            ],
        }
    )
    df_same_counts = assign_case_control_counts(df_same, case_samples, control_samples, all_samples)
    for col in [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]:
        assert col not in df_same_counts.columns

    res_same = perform_gene_burden_analysis(df_same_counts, cfg)
    assert len(res_same) == 1
    assert res_same.iloc[0]["proband_count"] == 1
    assert res_same.iloc[0]["control_count"] == 2

    # Case B: Different samples (P1 at site 1, P2 at site 2) missing -> 0 cases -> skipped
    df_diff = pd.DataFrame(
        {
            "GENE": ["GENE_DIFF", "GENE_DIFF"],
            "GT": [
                "P1(./.);P2(0/1);C1(0/0);C2(0/0)",
                "P1(0/0);P2(./.);C1(0/0);C2(0/0)",
            ],
        }
    )
    df_diff_counts = assign_case_control_counts(df_diff, case_samples, control_samples, all_samples)
    for col in [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]:
        assert col not in df_diff_counts.columns

    res_diff = perform_gene_burden_analysis(df_diff_counts, cfg)
    assert res_diff.empty


@pytest.mark.unit
def test_aggregate_gene_burden_preserves_cohort_proband_count():
    """Verify _aggregate_gene_burden_from_columns and _aggregate_gene_burden_from_gt
    preserve proband_count == len(case_samples) while exposing proband_called_count.
    """
    from variantcentrifuge.gene_burden import (
        _aggregate_gene_burden_from_columns,
        _aggregate_gene_burden_from_gt,
    )

    case_samples = {"Case1", "Case2"}
    control_samples = {"Ctrl1", "Ctrl2"}
    vcf_samples = ["Case1", "Case2", "Ctrl1", "Ctrl2"]
    gt_cols = ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT", "GEN_3__GT"]

    # 1 case missing (Case2 is ./.)
    df_cols = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "GEN_0__GT": ["0/1"],
            "GEN_1__GT": ["./."],
            "GEN_2__GT": ["0/0"],
            "GEN_3__GT": ["0/0"],
        }
    )
    agg_cols = _aggregate_gene_burden_from_columns(
        df_cols, case_samples, control_samples, vcf_samples, gt_cols
    )
    assert len(agg_cols) == 1
    assert agg_cols[0]["proband_count"] == 2  # Total cohort cases preserved for SKAT/regression
    assert agg_cols[0]["control_count"] == 2
    assert agg_cols[0]["proband_called_count"] == 1  # Called cases exposed for Fisher
    assert agg_cols[0]["control_called_count"] == 2

    df_gt = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "GT": ["Case1(0/1);Case2(./.);Ctrl1(0/0);Ctrl2(0/0)"],
        }
    )
    agg_gt = _aggregate_gene_burden_from_gt(df_gt, case_samples, control_samples)
    assert len(agg_gt) == 1
    assert agg_gt[0]["proband_count"] == 2
    assert agg_gt[0]["control_count"] == 2
    assert agg_gt[0]["proband_called_count"] == 1
    assert agg_gt[0]["control_called_count"] == 2


@pytest.mark.unit
def test_gt_aggregation_mixed_ploidy_order_independence():
    """Verify order independence when a sample has mixed ploidy calls across variants."""
    from variantcentrifuge.gene_burden import _aggregate_gene_burden_from_gt

    df_order1 = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "GT": ["Case1(1/1);Ctrl1(0/0)", "Case1(0);Ctrl1(0/0)"],
        }
    )
    df_order2 = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "GT": ["Case1(0);Ctrl1(0/0)", "Case1(1/1);Ctrl1(0/0)"],
        }
    )
    case_samples = {"Case1"}
    control_samples = {"Ctrl1"}
    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    agg1 = _aggregate_gene_burden_from_gt(df_order1, case_samples, control_samples)
    agg2 = _aggregate_gene_burden_from_gt(df_order2, case_samples, control_samples)
    assert agg1[0]["proband_called_chromosomes"] == 2
    assert agg2[0]["proband_called_chromosomes"] == 2
    assert agg1 == agg2

    res1 = perform_gene_burden_analysis(
        df_order1, cfg, case_samples=case_samples, control_samples=control_samples
    )
    res2 = perform_gene_burden_analysis(
        df_order2, cfg, case_samples=case_samples, control_samples=control_samples
    )

    assert len(res1) == 1
    assert len(res2) == 1
    assert res1.iloc[0]["raw_p_value"] == res2.iloc[0]["raw_p_value"]


@pytest.mark.unit
def test_sample_mixed_haploid_and_elided_diploid_reference():
    """Verify a sample with ['1', '0/0'] matches between raw-column mode
    and reconstructed GT mode.
    """
    from variantcentrifuge.gene_burden import (
        _aggregate_gene_burden_from_columns,
        _aggregate_gene_burden_from_gt,
    )

    # Raw column mode: Case1 has ['1', '0/0'], Ctrl1 has ['0/0', '0/0']
    df_cols = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "CHROM": ["chrX", "chrX"],
            "POS": [1000, 2000],
            "GEN_0__GT": ["1", "0/0"],
            "GEN_1__GT": ["0/0", "0/0"],
        }
    )
    case_samples = {"Case1"}
    control_samples = {"Ctrl1"}
    vcf_samples = ["Case1", "Ctrl1"]
    gt_cols = ["GEN_0__GT", "GEN_1__GT"]

    # Reconstructed GT mode: on row 2, '0/0' is elided
    df_gt = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "GT": ["Case1(1)", ""],
        }
    )

    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    res_cols = perform_gene_burden_analysis(
        df_cols,
        cfg,
        case_samples=case_samples,
        control_samples=control_samples,
        vcf_samples=vcf_samples,
    )
    res_gt = perform_gene_burden_analysis(
        df_gt, cfg, case_samples=case_samples, control_samples=control_samples
    )

    agg_cols = _aggregate_gene_burden_from_columns(
        df_cols, case_samples, control_samples, vcf_samples, gt_cols
    )
    agg_gt = _aggregate_gene_burden_from_gt(df_gt, case_samples, control_samples)

    assert agg_cols[0]["proband_called_chromosomes"] == 2
    assert agg_gt[0]["proband_called_chromosomes"] == 2
    assert res_cols.iloc[0]["raw_p_value"] == res_gt.iloc[0]["raw_p_value"]


@pytest.mark.unit
def test_legacy_aggregation_haploid_allele_cap():
    """Verify legacy aggregation caps allele count against proband_called_chromosomes."""
    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "proband_count": [1, 1],
            "control_count": [1, 1],
            "proband_called_chromosomes": [1, 1],
            "control_called_chromosomes": [2, 2],
            "proband_variant_count": [1, 1],
            "control_variant_count": [0, 0],
            "proband_allele_count": [1, 1],
            "control_allele_count": [0, 0],
        }
    )
    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}
    res = perform_gene_burden_analysis(df, cfg)
    assert len(res) == 1
    row = res.iloc[0]
    # 1 haploid case with 2 variants has sum(allele_count) = 2, capped at p_chrom = 1
    assert row["proband_allele_count"] == 1
    assert not np.isnan(row["raw_p_value"])


@pytest.mark.unit
def test_legacy_aggregation_excludes_incomplete_sample_from_numerators():
    """Verify that a sample missing at any site in a gene is excluded from carrier/allele counts
    in legacy aggregation (P1=[0/1, ./.], P2=[0/0, 0/0], C1=[0/1, 0/0], C2=[0/0, 0/0]).
    """
    from variantcentrifuge.gene_burden import (
        perform_gene_burden_analysis,
    )
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts

    case_samples = {"P1", "P2"}
    control_samples = {"C1", "C2"}
    all_samples = {"P1", "P2", "C1", "C2"}
    cfg = {"gene_burden_mode": "samples", "correction_method": "fdr"}

    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "GT": [
                "P1(0/1);P2(0/0);C1(0/1);C2(0/0)",
                "P1(./.);P2(0/0);C1(0/0);C2(0/0)",
            ],
        }
    )

    df_counts = assign_case_control_counts(df, case_samples, control_samples, all_samples)
    for col in [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]:
        assert col not in df_counts.columns

    # 3. perform_gene_burden_analysis with GT (sample-aware packed GT path)
    res_gt = perform_gene_burden_analysis(df_counts, cfg)
    assert len(res_gt) == 1
    assert res_gt.iloc[0]["proband_carrier_count"] == 0
    assert res_gt.iloc[0]["control_carrier_count"] == 1
    assert res_gt.iloc[0]["odds_ratio"] == pytest.approx(1.0 / 3.0)


@pytest.mark.unit
def test_haploid_multisite_chromosome_count_with_missingness():
    """Verify that multi-site haploid cases with incomplete samples compute
    exact complete-case called chromosome counts (proband_called_chromosomes == 1).
    """
    case_samples = {"P1", "P2", "P3"}
    control_samples = {"C1"}

    # P1 is (1, 1), P2 is (0, .), P3 is (., 0), C1 is (0/1, 0/1)
    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "CHROM": ["chrX", "chrX"],
            "POS": [1000, 2000],
            "GT": [
                "P1(1);P2(0);P3(.);C1(0/1)",
                "P1(1);P2(.);P3(0);C1(0/1)",
            ],
        }
    )

    # 1. Packed GT aggregation path
    gt_data = _aggregate_gene_burden_from_gt(df, case_samples, control_samples)
    assert len(gt_data) == 1
    assert gt_data[0]["proband_called_count"] == 1  # Only P1 is complete
    assert gt_data[0]["proband_called_chromosomes"] == 1  # P1 is haploid!
    assert gt_data[0]["control_called_count"] == 1
    assert gt_data[0]["control_called_chromosomes"] == 2  # C1 is diploid!


@pytest.mark.unit
def test_assign_case_control_counts_does_not_leak_tracking_columns():
    """Verify assign_case_control_counts does not leak internal sample tracking columns."""
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts

    case_samples = {"P_dip_ref", "P_hap_ref", "P_carrier"}
    control_samples = {"C_dip_ref", "C_carrier"}
    all_samples = case_samples | control_samples

    df = pd.DataFrame(
        {
            "GENE": ["GENE1"],
            "CHROM": ["chr1"],
            "POS": [1000],
            "GT": [
                "P_dip_ref(0/0);P_hap_ref(0);P_carrier(0/1);C_dip_ref(0/0);C_carrier(1/1)",
            ],
        }
    )
    df_counts = assign_case_control_counts(df, case_samples, control_samples, all_samples)
    for col in [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]:
        assert col not in df_counts.columns


@pytest.mark.unit
def test_colon_containing_sample_ids_in_gene_burden():
    """Verify that sample IDs containing colons (e.g. family:child) are parsed
    correctly in both packed-GT and legacy count paths.
    """
    from variantcentrifuge.gene_burden import perform_gene_burden_analysis
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts

    case_samples = {"family:child"}
    control_samples = {"family:parent"}
    all_samples = case_samples | control_samples
    cfg = {"gene_burden_mode": "alleles", "correction_method": "fdr"}

    df = pd.DataFrame(
        {
            "GENE": ["GENE1", "GENE1"],
            "CHROM": ["chr1", "chr1"],
            "POS": [1000, 2000],
            "GT": [
                "family:child(0/1);family:parent(0/0)",
                "family:child(0/0);family:parent(0/1)",
            ],
        }
    )
    df_counts = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # 1. Packed GT path (alleles mode)
    res_gt = perform_gene_burden_analysis(df_counts, cfg)
    assert len(res_gt) == 1
    assert res_gt.iloc[0]["proband_allele_count"] == 1
    assert res_gt.iloc[0]["control_allele_count"] == 1

    # 2. Legacy / counts-only path (alleles mode)
    df_no_gt = df_counts.drop(columns=["GT"])
    res_legacy = perform_gene_burden_analysis(df_no_gt, cfg)
    assert len(res_legacy) == 1
    assert res_legacy.iloc[0]["proband_allele_count"] == 1
    assert res_legacy.iloc[0]["control_allele_count"] == 1

    # 3. Samples mode (carrier counts)
    cfg_samp = {"gene_burden_mode": "samples", "correction_method": "fdr"}
    res_gt_samp = perform_gene_burden_analysis(df_counts, cfg_samp)
    assert res_gt_samp.iloc[0]["proband_carrier_count"] == 1
    assert res_gt_samp.iloc[0]["control_carrier_count"] == 1
    res_leg_samp = perform_gene_burden_analysis(df_no_gt, cfg_samp)
    assert res_leg_samp.iloc[0]["proband_variant_count"] == 1
    assert res_leg_samp.iloc[0]["control_variant_count"] == 1


@pytest.mark.unit
def test_pseudonymization_and_variant_level_output_tracking_columns():
    """Verify that internal tracking columns are omitted from variant-level TSV output
    and that SamplePseudonymizer never leaks original sample IDs.
    """
    from variantcentrifuge.analyze_variants import analyze_variants
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts
    from variantcentrifuge.pseudonymizer import create_pseudonymizer

    case_samples = {"PATIENT_A"}
    control_samples = {"PATIENT_B"}
    all_samples = case_samples | control_samples

    df = pd.DataFrame(
        {
            "CHROM": ["chr1"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["G"],
            "GENE": ["GENE1"],
            "FILTER": ["PASS"],
            "GT": ["PATIENT_A(0/1);PATIENT_B(.)"],
        }
    )

    df_counts = assign_case_control_counts(df, case_samples, control_samples, all_samples)
    tracking_columns = [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]
    for col in tracking_columns:
        assert col not in df_counts.columns

    # 1. Verify pseudonymize_dataframe pseudonymizes internal tracking columns if present
    pseudonymizer = create_pseudonymizer("sequential", prefix="SAMPLE")
    pseudonymizer.create_mapping(sorted(all_samples))
    df_with_tracking = df_counts.copy()
    df_with_tracking["proband_carrier_samples"] = ["PATIENT_A:1:2"]
    df_with_tracking["control_missing_samples"] = ["PATIENT_B"]
    pseudo_df = pseudonymizer.pseudonymize_dataframe(df_with_tracking)

    assert "PATIENT_A" not in pseudo_df["proband_carrier_samples"].iloc[0]
    assert "PATIENT_B" not in pseudo_df["control_missing_samples"].iloc[0]
    assert "SAMPLE_001" in pseudo_df["proband_carrier_samples"].iloc[0]
    assert "SAMPLE_002" in pseudo_df["control_missing_samples"].iloc[0]

    # 2. Verify analyze_variants omits tracking columns from variant output lines
    tsv_lines = df_counts.to_csv(sep="\t", index=False).splitlines(keepends=True)
    cfg = {
        "sample_list": "PATIENT_A,PATIENT_B",
        "perform_gene_burden": False,
        "no_stats": True,
    }
    lines = list(analyze_variants(tsv_lines, cfg))
    assert len(lines) == 2

    tracking_columns = [
        "proband_missing_samples",
        "control_missing_samples",
        "proband_carrier_samples",
        "control_carrier_samples",
    ]
    header = lines[0].split("\t")
    for col in tracking_columns:
        assert col not in header

    variant_line = lines[1].split("\t")
    assert len(variant_line) == len(header)
    assert "chr1" in variant_line
    assert "GENE1" in variant_line
    for col in tracking_columns:
        assert col not in variant_line


def test_assign_case_control_counts_exact_called_denominators():
    """Absent samples count as called 0/0; listed missing calls are not called (bulk path)."""
    from variantcentrifuge.phenotype_helpers import assign_case_control_counts

    # E is in both sets and counts as a case; X is not in the VCF and is ignored.
    df = pd.DataFrame({"GT": ["A(./.);B(1/1);C(1);X(0/1)", ""]})
    out = assign_case_control_counts(
        df, {"A", "B", "E"}, {"C", "D", "E"}, {"A", "B", "C", "D", "E"}
    )
    cols = [
        "proband_called_count",
        "proband_called_chromosomes",
        "proband_variant_count",
        "proband_allele_count",
        "proband_homozygous_count",
        "control_called_count",
        "control_called_chromosomes",
        "control_variant_count",
        "control_allele_count",
        "control_homozygous_count",
    ]
    # Row 0: cases B (1/1) + E (implicit 0/0), A missing; controls C (haploid 1) + D (0/0).
    assert out.loc[0, cols].tolist() == [2, 4, 1, 2, 1, 2, 3, 1, 1, 0]
    # Row 1: nobody listed, so every sample is a called 0/0.
    assert out.loc[1, cols].tolist() == [3, 6, 0, 0, 0, 2, 4, 0, 0, 0]

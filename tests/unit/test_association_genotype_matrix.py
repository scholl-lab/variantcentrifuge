"""
Unit tests for variantcentrifuge.association.genotype_matrix.

Tests cover:
- parse_gt_to_dosage: standard, missing, phased, multi-allelic, partial missing (BURDEN-03)
- build_genotype_matrix: shape, imputation, site filter, sample filter, MAF, warnings
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.genotype_matrix import build_genotype_matrix, parse_gt_to_dosage

# ---------------------------------------------------------------------------
# parse_gt_to_dosage tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestParseGtToDosage:
    """parse_gt_to_dosage edge cases — BURDEN-03."""

    # Standard unphased genotypes

    def test_parse_gt_hom_ref(self) -> None:
        """0/0 -> dosage=0, is_multi_allelic=False."""
        dosage, is_multi = parse_gt_to_dosage("0/0")
        assert dosage == 0
        assert is_multi is False

    def test_parse_gt_het_01(self) -> None:
        """0/1 -> dosage=1, is_multi_allelic=False."""
        dosage, is_multi = parse_gt_to_dosage("0/1")
        assert dosage == 1
        assert is_multi is False

    def test_parse_gt_het_10(self) -> None:
        """1/0 -> dosage=1, is_multi_allelic=False (order doesn't matter)."""
        dosage, is_multi = parse_gt_to_dosage("1/0")
        assert dosage == 1
        assert is_multi is False

    def test_parse_gt_hom_alt(self) -> None:
        """1/1 -> dosage=2, is_multi_allelic=False."""
        dosage, is_multi = parse_gt_to_dosage("1/1")
        assert dosage == 2
        assert is_multi is False

    # Missing genotypes

    def test_parse_gt_missing_dot_slash_dot(self) -> None:
        """./. -> (None, False)."""
        dosage, is_multi = parse_gt_to_dosage("./.")
        assert dosage is None
        assert is_multi is False

    def test_parse_gt_missing_phased(self) -> None:
        """.|. -> (None, False)."""
        dosage, is_multi = parse_gt_to_dosage(".|.")
        assert dosage is None
        assert is_multi is False

    def test_parse_gt_missing_empty_string(self) -> None:
        """Empty string -> (None, False)."""
        dosage, is_multi = parse_gt_to_dosage("")
        assert dosage is None
        assert is_multi is False

    def test_parse_gt_missing_none_input(self) -> None:
        """None input -> (None, False)."""
        dosage, is_multi = parse_gt_to_dosage(None)
        assert dosage is None
        assert is_multi is False

    # Phased genotypes

    def test_parse_gt_phased_hom_ref(self) -> None:
        """0|0 -> (0, False)."""
        dosage, is_multi = parse_gt_to_dosage("0|0")
        assert dosage == 0
        assert is_multi is False

    def test_parse_gt_phased_het_01(self) -> None:
        """0|1 -> (1, False)."""
        dosage, is_multi = parse_gt_to_dosage("0|1")
        assert dosage == 1
        assert is_multi is False

    def test_parse_gt_phased_het_10(self) -> None:
        """1|0 -> (1, False)."""
        dosage, is_multi = parse_gt_to_dosage("1|0")
        assert dosage == 1
        assert is_multi is False

    def test_parse_gt_phased_hom_alt(self) -> None:
        """1|1 -> (2, False)."""
        dosage, is_multi = parse_gt_to_dosage("1|1")
        assert dosage == 2
        assert is_multi is False

    # Multi-allelic genotypes — critical BURDEN-03 test

    def test_parse_gt_multi_allelic_1_2(self) -> None:
        """1/2 -> dosage=1 (het-equivalent), is_multi_allelic=True."""
        dosage, is_multi = parse_gt_to_dosage("1/2")
        assert dosage == 1
        assert is_multi is True

    def test_parse_gt_multi_allelic_0_2(self) -> None:
        """0/2 -> dosage=1, is_multi_allelic=True."""
        dosage, is_multi = parse_gt_to_dosage("0/2")
        assert dosage == 1
        assert is_multi is True

    def test_parse_gt_multi_allelic_2_2(self) -> None:
        """2/2 -> dosage=1, is_multi_allelic=True."""
        dosage, is_multi = parse_gt_to_dosage("2/2")
        assert dosage == 1
        assert is_multi is True

    def test_parse_gt_multi_allelic_phased(self) -> None:
        """1|2 -> (1, True) — phased multi-allelic."""
        dosage, is_multi = parse_gt_to_dosage("1|2")
        assert dosage == 1
        assert is_multi is True

    # Partial missing

    def test_parse_gt_partial_missing_dot_one(self) -> None:
        """./1 -> (None, False) — partial missing is treated as missing."""
        dosage, is_multi = parse_gt_to_dosage("./1")
        assert dosage is None
        assert is_multi is False

    def test_parse_gt_partial_missing_one_dot(self) -> None:
        """1/. -> (None, False) — partial missing is treated as missing."""
        dosage, is_multi = parse_gt_to_dosage("1/.")
        assert dosage is None
        assert is_multi is False

    def test_parse_gt_all_standard_unphased(self) -> None:
        """Batch check: all standard unphased GTs produce correct (dosage, False)."""
        expected = {
            "0/0": (0, False),
            "0/1": (1, False),
            "1/0": (1, False),
            "1/1": (2, False),
        }
        for gt, (exp_dosage, exp_multi) in expected.items():
            dosage, is_multi = parse_gt_to_dosage(gt)
            assert dosage == exp_dosage, f"GT={gt}: expected dosage {exp_dosage}, got {dosage}"
            assert is_multi is exp_multi, f"GT={gt}: expected is_multi={exp_multi}, got {is_multi}"


# ---------------------------------------------------------------------------
# build_genotype_matrix helper
# ---------------------------------------------------------------------------


def _make_gene_df(
    n_samples: int,
    n_variants: int,
    gt_values: list[list[str]] | None = None,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """
    Build a minimal gene DataFrame for testing build_genotype_matrix.

    Parameters
    ----------
    n_samples : int
        Number of samples.
    n_variants : int
        Number of variants (rows).
    gt_values : list of lists, optional
        Outer list is per-variant, inner per-sample. Default: all "0/0".

    Returns
    -------
    gene_df, vcf_samples, gt_columns
    """
    vcf_samples = [f"SAMPLE_{i}" for i in range(n_samples)]
    gt_columns = [f"GEN_{i}__GT" for i in range(n_samples)]

    if gt_values is None:
        gt_values = [["0/0"] * n_samples for _ in range(n_variants)]

    data = {}
    for s_idx, col in enumerate(gt_columns):
        data[col] = [gt_values[v][s_idx] for v in range(n_variants)]

    gene_df = pd.DataFrame(data)
    return gene_df, vcf_samples, gt_columns


# ---------------------------------------------------------------------------
# build_genotype_matrix tests
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestBuildGenotypeMatrixBasic:
    """Basic shape, dtype, and output structure tests for build_genotype_matrix."""

    def test_build_genotype_matrix_basic_shape(self) -> None:
        """3 samples, 2 variants, no missing -> shape (3, 2), no NaN."""
        gt_values = [
            ["1/1", "0/0", "0/1"],
            ["0/0", "0/1", "0/0"],
        ]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 2, gt_values)

        geno, mafs, sample_mask, _ = build_genotype_matrix(
            gene_df, vcf_samples, gt_cols, is_binary=True
        )

        assert geno.shape == (3, 2), f"Expected (3, 2), got {geno.shape}"
        assert geno.dtype == np.float64
        assert not np.isnan(geno).any(), "No NaN expected with no missing data"
        assert len(sample_mask) == 3
        assert len(mafs) == 2

    def test_build_genotype_matrix_dosage_values(self) -> None:
        """Dosages are correctly mapped: 0/0->0, 0/1->1, 1/1->2."""
        gt_values = [
            ["0/0", "0/1", "1/1"],
        ]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 1, gt_values)

        geno, _, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols, is_binary=True)

        assert geno[0, 0] == pytest.approx(0.0)  # 0/0
        assert geno[1, 0] == pytest.approx(1.0)  # 0/1
        assert geno[2, 0] == pytest.approx(2.0)  # 1/1

    def test_build_genotype_matrix_maf_computation(self) -> None:
        """MAF computed correctly from observed dosages."""
        # 4 samples: 0/0, 0/1, 0/1, 1/1
        # Dosages: 0, 1, 1, 2 -> mean = 4/4 = 1.0 -> MAF = 0.5
        gt_values = [["0/0", "0/1", "0/1", "1/1"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(4, 1, gt_values)

        _, mafs, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert mafs[0] == pytest.approx(0.5), f"Expected MAF=0.5, got {mafs[0]}"

    def test_build_genotype_matrix_maf_rare_variant(self) -> None:
        """Rare variant (1 het in 10 samples) has MAF=0.05."""
        gt_values = [["0/1"] + ["0/0"] * 9]
        gene_df, vcf_samples, gt_cols = _make_gene_df(10, 1, gt_values)

        _, mafs, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        # mean dosage = 1/10 = 0.1 -> MAF = 0.1/2 = 0.05
        assert mafs[0] == pytest.approx(0.05), f"Expected MAF=0.05, got {mafs[0]}"

    def test_build_genotype_matrix_sample_mask_all_true_no_missing(self) -> None:
        """sample_mask is all True when no samples have missing data."""
        gt_values = [["0/1", "0/0", "1/1"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 1, gt_values)

        _, _, sample_mask, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert all(sample_mask), f"Expected all True, got {sample_mask}"


@pytest.mark.unit
class TestBuildGenotypeMatrixMissingData:
    """Missing genotype handling: site filter, sample filter, imputation."""

    def test_build_genotype_matrix_missing_imputation_binary(self) -> None:
        """Binary trait: missing GTs are imputed with round(2*MAF)."""
        # 20 samples: S0=0/1, S1=./., S2-S19=0/0
        # Missing rate = 1/20 = 5% (below 10% threshold -> variant kept)
        n_samples = 20
        gt_values = [["0/1", "./."] + ["0/0"] * (n_samples - 2)]
        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, 1, gt_values)

        geno, _, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols, is_binary=True)

        # Variant should be kept (5% missing < 10% threshold)
        assert geno.shape[1] == 1, "Variant should be kept (5% missing < threshold)"
        # No NaN after imputation
        assert not np.isnan(geno).any()
        # Imputed value at S1 (index 1) should be an integer (0, 1, or 2)
        imputed_val = geno[1, 0]
        assert imputed_val in (0.0, 1.0, 2.0), (
            f"Binary imputed value should be 0/1/2, got {imputed_val}"
        )

    def test_build_genotype_matrix_missing_imputation_continuous(self) -> None:
        """Continuous trait: missing GTs are imputed with 2*MAF (may be fractional)."""
        # 20 samples: S0=0/1, S1=./., S2..S10=0/0 (9), S11..S19=0/1 (9)
        # Missing rate = 1/20 = 5% (below 10% threshold -> variant kept)
        n_samples = 20
        gt_values = [["0/1", "./."] + ["0/0"] * 9 + ["0/1"] * 9]
        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, 1, gt_values)

        geno, mafs, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols, is_binary=False)

        assert geno.shape[1] == 1, "Variant should be kept (5% missing < threshold)"
        assert not np.isnan(geno).any()
        imputed_val = geno[1, 0]
        expected = 2.0 * mafs[0]
        assert imputed_val == pytest.approx(expected), (
            f"Continuous imputed value should be 2*MAF={expected}, got {imputed_val}"
        )

    def test_build_genotype_matrix_site_filter_removes_high_missing_variant(self) -> None:
        """Variant with >10% missing is removed from output matrix."""
        n_samples = 20
        # Variant 0: all valid (0/0); Variant 1: 50% missing (>10% -> removed)
        gt_variant_0 = ["0/0"] * n_samples
        gt_variant_1 = ["./."] * (n_samples // 2) + ["0/0"] * (n_samples // 2)

        gt_values = [gt_variant_0, gt_variant_1]
        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, 2, gt_values)

        geno, mafs, _, _ = build_genotype_matrix(
            gene_df, vcf_samples, gt_cols, missing_site_threshold=0.10
        )

        # Only 1 variant should survive (variant 1 removed)
        assert geno.shape[1] == 1, f"Expected 1 kept variant, got {geno.shape[1]}"
        assert len(mafs) == 1

    def test_build_genotype_matrix_site_filter_keeps_within_threshold(self) -> None:
        """Variant with exactly 10% missing (at threshold) is kept."""
        n_samples = 10
        # 1 of 10 missing = 10% -> at threshold (<=10% kept)
        gt_values = [["./."] + ["0/1"] * 9]
        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, 1, gt_values)

        geno, _, _, _ = build_genotype_matrix(
            gene_df, vcf_samples, gt_cols, missing_site_threshold=0.10
        )

        # Variant should be kept (1/10 = 10% <= threshold)
        assert geno.shape[1] == 1, "Variant at exactly threshold should be kept"

    def test_build_genotype_matrix_sample_filter_flags_high_missing(self) -> None:
        """Sample with >80% missing across variants is flagged in sample_mask."""
        # S0 is missing for all variants (100% missing per sample).
        # Each variant has 1/5=20% missing -> use threshold=0.25 to keep variants.
        # Then sample filter: S0 has 100% missing > 80% threshold -> flagged False.
        n_variants = 10
        n_samples = 5

        gt_values = []
        for _ in range(n_variants):
            row = ["./."] + ["0/1"] * (n_samples - 1)
            gt_values.append(row)

        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, n_variants, gt_values)

        _, _, sample_mask, _ = build_genotype_matrix(
            gene_df,
            vcf_samples,
            gt_cols,
            missing_site_threshold=0.25,  # 20% missing per site < 25% -> kept
            missing_sample_threshold=0.80,
        )

        assert len(sample_mask) == n_samples
        # S0 has 100% missing -> should be flagged as False
        assert sample_mask[0] is False, "High-missing sample should have False in sample_mask"
        # Other samples should be True
        assert all(sample_mask[1:]), "Low-missing samples should have True in sample_mask"

    def test_build_genotype_matrix_no_nan_after_imputation(self) -> None:
        """Returned geno matrix has zero NaN entries after imputation."""
        # Mix of missing, standard, phased genotypes
        gt_values = [
            ["0/1", "./.", "1/1", "0|1", "0/0"],
            ["./.", "0/0", "0/1", "./.", "1/1"],
        ]
        gene_df, vcf_samples, gt_cols = _make_gene_df(5, 2, gt_values)

        geno, _, _, _ = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert not np.isnan(geno).any(), "geno must contain no NaN after build_genotype_matrix"


@pytest.mark.unit
class TestBuildGenotypeMatrixMAF:
    """MAF computation tests — must use all samples, not phenotype-stratified."""

    def test_build_genotype_matrix_maf_from_all_samples(self) -> None:
        """MAF is computed from all samples combined, not stratified by phenotype."""
        # 6 samples: 3 cases (S0-S2), 3 controls (S3-S5)
        # Variant: S0=0/1, S1=0/0, S2=0/0, S3=0/0, S4=0/0, S5=0/0
        # All-sample MAF = (1.0 / 6.0) / 2.0 = 0.0833
        # Case-only MAF would be (1 / 3) / 2 = 0.167 (different from all-sample)
        gt_values = [["0/1", "0/0", "0/0", "0/0", "0/0", "0/0"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(6, 1, gt_values)

        phenotype = np.array([1, 1, 1, 0, 0, 0], dtype=float)  # cases=S0-S2

        _, mafs, _, _ = build_genotype_matrix(
            gene_df, vcf_samples, gt_cols, phenotype_vector=phenotype
        )

        # MAF from all samples: mean dosage = 1/6, MAF = mean/2 = 1/12
        expected_maf = (1.0 / 6.0) / 2.0
        assert mafs[0] == pytest.approx(expected_maf, abs=1e-10), (
            f"MAF should be computed from all samples: expected {expected_maf}, got {mafs[0]}"
        )

    def test_build_genotype_matrix_every_sample_missing_maf_zero(self) -> None:
        """All-missing variant gets MAF=0 and imputed to 0."""
        gt_values = [["./.", "./.", "./."]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 1, gt_values)

        # Variant has 100% missing — use high threshold to keep it
        geno, mafs, _, _ = build_genotype_matrix(
            gene_df, vcf_samples, gt_cols, missing_site_threshold=1.0
        )

        assert mafs[0] == 0.0
        assert (geno[:, 0] == 0.0).all()


@pytest.mark.unit
class TestBuildGenotypeMatrixWarnings:
    """Warning message generation tests."""

    def test_build_genotype_matrix_multi_allelic_warning(self) -> None:
        """1/2 genotype triggers multi-allelic warning in warnings_list."""
        gt_values = [["1/2", "0/0", "0/1"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 1, gt_values)

        _, _, _, warnings_list = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert len(warnings_list) > 0, "Expected multi-allelic warning"
        assert any(
            "multi-allelic" in w.lower() or "bcftools" in w.lower() for w in warnings_list
        ), f"Expected bcftools norm recommendation in warnings, got: {warnings_list}"

    def test_build_genotype_matrix_no_warnings_clean_data(self) -> None:
        """Clean data (no multi-allelic, no differential missingness) -> empty warnings."""
        gt_values = [["0/0", "0/1", "1/1"], ["0/1", "0/0", "0/0"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 2, gt_values)

        _, _, _, warnings_list = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert warnings_list == [], f"Expected no warnings, got: {warnings_list}"

    def test_build_genotype_matrix_multiple_multi_allelic_counted(self) -> None:
        """Multiple multi-allelic genotypes are counted in the warning message."""
        # 2 variants, each has a 1/2 genotype -> 2 total multi-allelic genotypes
        gt_values = [
            ["1/2", "0/0", "0/1"],
            ["0/0", "1/2", "0/1"],
        ]
        gene_df, vcf_samples, gt_cols = _make_gene_df(3, 2, gt_values)

        _, _, _, warnings_list = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert len(warnings_list) > 0
        # The warning should mention the count (2 genotypes)
        warning_text = " ".join(warnings_list)
        assert "2" in warning_text, f"Expected count '2' in warning, got: {warning_text}"

    def test_build_genotype_matrix_returns_four_tuple(self) -> None:
        """Return value is a 4-tuple (geno, mafs, sample_mask, warnings_list)."""
        gt_values = [["0/1", "0/0"]]
        gene_df, vcf_samples, gt_cols = _make_gene_df(2, 1, gt_values)

        result = build_genotype_matrix(gene_df, vcf_samples, gt_cols)

        assert isinstance(result, tuple)
        assert len(result) == 4
        geno, mafs, sample_mask, warnings_list = result
        assert isinstance(geno, np.ndarray)
        assert isinstance(mafs, np.ndarray)
        assert isinstance(sample_mask, list)
        assert isinstance(warnings_list, list)

    def test_build_genotype_matrix_differential_missingness_warning(self) -> None:
        """Differential missingness >5% between cases/controls triggers warning."""
        # 10 samples: 5 cases, 5 controls
        # Variant 0: all cases missing (100%), controls have data (0% missing)
        n_samples = 10
        gt_values = [
            ["./."] * 5 + ["0/1"] * 5,  # cases all missing
        ]
        gene_df, vcf_samples, gt_cols = _make_gene_df(n_samples, 1, gt_values)

        phenotype = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0], dtype=float)

        _, _, _, warnings_list = build_genotype_matrix(
            gene_df,
            vcf_samples,
            gt_cols,
            missing_site_threshold=1.0,  # high threshold so variant is kept
            phenotype_vector=phenotype,
        )

        # At least one warning about differential missingness
        assert any("differential" in w.lower() for w in warnings_list), (
            f"Expected differential missingness warning, got: {warnings_list}"
        )

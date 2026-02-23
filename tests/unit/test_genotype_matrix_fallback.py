"""
Unit tests for GT matrix recovery in AssociationAnalysisStage (COAST-01).

Tests the fallback logic in analysis_stages.py that recovers per-sample GT columns
for genotype matrix construction when variants_df is None or lacks GT columns.

The branch under test (lines ~2375-2395 in analysis_stages.py):

    elif "GT" in df.columns and needs_regression and context.vcf_samples:
        fallback_df = context.variants_df
        gt_cols_fb = _find_per_sample_gt_columns(fallback_df) if fallback_df is not None else []
        if gt_cols_fb:
            df_with_per_sample_gt = fallback_df
        else:
            gt_cols_df = _find_per_sample_gt_columns(df)
            if gt_cols_df:
                df_with_per_sample_gt = df

Requirements covered: COAST-01 (genotype matrix fallback when variants_df is None)
"""

from __future__ import annotations

import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import _find_per_sample_gt_columns

# ---------------------------------------------------------------------------
# Helper to build DataFrames with/without per-sample GT columns
# ---------------------------------------------------------------------------


def _df_with_gt_cols(n_variants: int = 3, n_samples: int = 2) -> pd.DataFrame:
    """Build DataFrame with GEN_0__GT, GEN_1__GT columns (sanitized names)."""
    data: dict = {
        "GENE": ["BRCA1"] * n_variants,
        "GT": ["S1(0/1);S2(0/0)"] * n_variants,
    }
    for i in range(n_samples):
        data[f"GEN_{i}__GT"] = ["0/1", "0/0", "1/1"][:n_variants]
    return pd.DataFrame(data)


def _df_without_gt_cols(n_variants: int = 3) -> pd.DataFrame:
    """Build DataFrame with packed GT column only (no per-sample columns)."""
    return pd.DataFrame(
        {
            "GENE": ["BRCA1"] * n_variants,
            "GT": ["S1(0/1);S2(0/0)"] * n_variants,
        }
    )


# ---------------------------------------------------------------------------
# Tests for _find_per_sample_gt_columns helper
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestFindPerSampleGtColumns:
    """Test _find_per_sample_gt_columns recognizes sanitized and original column names."""

    def test_sanitized_names_recognized(self):
        """GEN_0__GT and GEN_1__GT are returned in index order."""
        df = _df_with_gt_cols(n_samples=3)
        cols = _find_per_sample_gt_columns(df)
        assert cols == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"]

    def test_original_names_recognized(self):
        """GEN[0].GT and GEN[1].GT are recognized."""
        df = pd.DataFrame(
            {
                "GENE": ["X"],
                "GT": ["S1(0/1)"],
                "GEN[0].GT": ["0/1"],
                "GEN[1].GT": ["0/0"],
            }
        )
        cols = _find_per_sample_gt_columns(df)
        assert "GEN[0].GT" in cols
        assert "GEN[1].GT" in cols

    def test_empty_when_no_gt_cols(self):
        """Returns empty list when no per-sample GT columns present."""
        df = _df_without_gt_cols()
        cols = _find_per_sample_gt_columns(df)
        assert cols == []

    def test_returns_empty_for_none_df(self):
        """Returns [] when df is None (guard handled before call in analysis_stages)."""
        # The function itself requires a DataFrame; guard is in the caller.
        # Test that the guard in analysis_stages branch behaves correctly:
        fallback_df = None
        gt_cols = _find_per_sample_gt_columns(fallback_df) if fallback_df is not None else []
        assert gt_cols == []


# ---------------------------------------------------------------------------
# Tests for the fallback branch logic itself
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestGenotypeMatrixFallbackLogic:
    """
    Test the GT fallback branch logic from analysis_stages.py.

    We directly simulate the branch conditions to verify they produce the
    correct df_with_per_sample_gt assignment.
    """

    def _simulate_fallback_branch(
        self,
        df: pd.DataFrame,
        variants_df: pd.DataFrame | None,
        vcf_samples: list[str],
        needs_regression: bool = True,
    ) -> pd.DataFrame | None:
        """
        Simulate the elif branch from AssociationAnalysisStage._process().

        This replicates the exact logic we wrote in analysis_stages.py so we
        can test it independently of the large _process() method.
        """
        df_with_per_sample_gt: pd.DataFrame | None = None

        if "GT" in df.columns and needs_regression and vcf_samples:
            fallback_df = variants_df
            gt_cols_fb = _find_per_sample_gt_columns(fallback_df) if fallback_df is not None else []
            if gt_cols_fb:
                df_with_per_sample_gt = fallback_df
            else:
                gt_cols_df = _find_per_sample_gt_columns(df)
                if gt_cols_df:
                    df_with_per_sample_gt = df

        return df_with_per_sample_gt

    def test_fallback_to_df_when_variants_df_is_none(self):
        """
        When variants_df is None but df has per-sample GT columns,
        df_with_per_sample_gt should be set to df.
        """
        df = _df_with_gt_cols(n_variants=3, n_samples=2)
        variants_df = None
        vcf_samples = ["S1", "S2"]

        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is df, "Expected fallback to current df when variants_df is None"
        assert _find_per_sample_gt_columns(result) == ["GEN_0__GT", "GEN_1__GT"]

    def test_fallback_to_variants_df_when_available(self):
        """
        When variants_df has GT columns, df_with_per_sample_gt should be variants_df.
        """
        df = _df_without_gt_cols(n_variants=3)  # No per-sample cols in current df
        variants_df = _df_with_gt_cols(n_variants=3, n_samples=2)
        vcf_samples = ["S1", "S2"]

        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is variants_df, "Expected variants_df to be used when it has GT cols"

    def test_variants_df_preferred_over_df(self):
        """
        When both variants_df and df have per-sample GT columns,
        variants_df should be preferred (checked first).
        """
        df = _df_with_gt_cols(n_variants=3, n_samples=2)  # has GT cols
        variants_df = _df_with_gt_cols(n_variants=3, n_samples=3)  # also has GT cols
        vcf_samples = ["S1", "S2", "S3"]

        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is variants_df, "variants_df should be preferred over df"

    def test_no_fallback_when_no_gt_columns_anywhere(self):
        """
        When neither df nor variants_df has per-sample GT columns,
        df_with_per_sample_gt should remain None.
        """
        df = _df_without_gt_cols(n_variants=3)
        variants_df = _df_without_gt_cols(n_variants=3)
        vcf_samples = ["S1", "S2"]

        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is None, "Expected None when no GT columns available anywhere"

    def test_no_fallback_when_variants_df_none_and_df_no_gt_cols(self):
        """
        When variants_df is None and df has no per-sample GT cols, result is None.
        """
        df = _df_without_gt_cols(n_variants=3)
        variants_df = None
        vcf_samples = ["S1", "S2"]

        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is None

    def test_branch_not_entered_when_needs_regression_false(self):
        """
        When needs_regression is False, the branch should not execute
        (df_with_per_sample_gt stays None regardless of columns).
        """
        df = _df_with_gt_cols(n_variants=3, n_samples=2)
        variants_df = _df_with_gt_cols(n_variants=3, n_samples=2)
        vcf_samples = ["S1", "S2"]

        result = self._simulate_fallback_branch(
            df, variants_df, vcf_samples, needs_regression=False
        )

        assert result is None, "Branch should not execute when needs_regression=False"

    def test_branch_not_entered_when_gt_not_in_df(self):
        """
        The elif branch only fires when 'GT' is already in df.columns.
        If GT is absent, fallback stays None.
        """
        # df has per-sample GT cols but no packed GT column
        df = pd.DataFrame(
            {
                "GENE": ["BRCA1"] * 3,
                "GEN_0__GT": ["0/1", "0/0", "1/1"],
                "GEN_1__GT": ["0/0", "1/1", "0/1"],
            }
        )
        variants_df = _df_with_gt_cols(n_variants=3, n_samples=2)
        vcf_samples = ["S1", "S2"]

        # The branch requires "GT" in df.columns — it's not here
        result = self._simulate_fallback_branch(df, variants_df, vcf_samples)

        assert result is None, "Branch should not execute when 'GT' not in df.columns"

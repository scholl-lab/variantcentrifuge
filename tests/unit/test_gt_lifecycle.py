"""Tests for GT column lifecycle through analysis stages (Fix 5)."""
import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    _find_per_sample_gt_columns,
    reconstruct_gt_column,
)


@pytest.mark.unit
class TestFindPerSampleGtColumns:
    """Verify the canonical GT column finder works for both name formats."""

    def test_sanitized_names(self):
        df = pd.DataFrame({"GEN_0__GT": ["0/1"], "GEN_1__GT": ["1/1"], "GENE": ["A"]})
        assert _find_per_sample_gt_columns(df) == ["GEN_0__GT", "GEN_1__GT"]

    def test_original_names(self):
        df = pd.DataFrame({"GEN[0].GT": ["0/1"], "GEN[1].GT": ["1/1"], "GENE": ["A"]})
        assert _find_per_sample_gt_columns(df) == ["GEN[0].GT", "GEN[1].GT"]

    def test_no_gt_columns(self):
        df = pd.DataFrame({"GENE": ["A"], "CHROM": ["1"]})
        assert _find_per_sample_gt_columns(df) == []

    def test_sorted_by_index(self):
        df = pd.DataFrame({"GEN_10__GT": ["0/1"], "GEN_2__GT": ["1/1"], "GEN_0__GT": ["0/0"]})
        assert _find_per_sample_gt_columns(df) == ["GEN_0__GT", "GEN_2__GT", "GEN_10__GT"]

    def test_mixed_with_non_gt_columns(self):
        df = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [100],
                "GEN_0__GT": ["0/1"],
                "GEN_1__GT": ["1/1"],
                "GENE": ["A"],
                "GEN_0__DP": [30],  # depth column — should NOT be matched
            }
        )
        result = _find_per_sample_gt_columns(df)
        assert result == ["GEN_0__GT", "GEN_1__GT"]
        assert "GEN_0__DP" not in result


@pytest.mark.unit
class TestReconstructGtColumnLocalCopy:
    """Verify that reconstruct_gt_column on a copy does not affect the original."""

    def test_original_preserved(self):
        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/1", "1/1"],
                "GEN_1__GT": ["0/0", "0/1"],
                "GENE": ["A", "B"],
            }
        )
        original_cols = set(df.columns)

        # Reconstruct on a copy (the Fix 5 pattern)
        reconstructed = reconstruct_gt_column(df.copy(), ["Sample1", "Sample2"])

        # Original unchanged
        assert set(df.columns) == original_cols
        assert "GEN_0__GT" in df.columns
        assert "GEN_1__GT" in df.columns

        # Reconstructed has packed GT, no per-sample cols
        assert "GT" in reconstructed.columns
        assert "GEN_0__GT" not in reconstructed.columns

    def test_reconstructed_has_packed_format(self):
        """Verify the packed GT format contains sample names and genotypes."""
        df = pd.DataFrame(
            {
                "GEN_0__GT": ["0/1"],
                "GEN_1__GT": ["1/1"],
                "GENE": ["A"],
            }
        )
        reconstructed = reconstruct_gt_column(df.copy(), ["CaseA", "ControlB"])
        assert "GT" in reconstructed.columns
        gt_val = reconstructed["GT"].iloc[0]
        assert "CaseA" in gt_val
        assert "ControlB" in gt_val


@pytest.mark.unit
class TestGtColumnsConsolidation:
    """Verify _find_gt_columns from gene_burden uses the same canonical implementation."""

    def test_gene_burden_find_gt_is_same_function(self):
        from variantcentrifuge.gene_burden import _find_gt_columns

        # After consolidation, _find_gt_columns should delegate to _find_per_sample_gt_columns
        df = pd.DataFrame({"GEN_0__GT": ["0/1"], "GEN_1__GT": ["1/1"]})
        assert _find_gt_columns(df) == _find_per_sample_gt_columns(df)

    def test_gene_burden_find_gt_sanitized(self):
        from variantcentrifuge.gene_burden import _find_gt_columns

        df = pd.DataFrame(
            {"GEN_0__GT": ["0/1"], "GEN_1__GT": ["0/0"], "GEN_2__GT": ["1/1"], "GENE": ["A"]}
        )
        assert _find_gt_columns(df) == ["GEN_0__GT", "GEN_1__GT", "GEN_2__GT"]

    def test_gene_burden_find_gt_empty(self):
        from variantcentrifuge.gene_burden import _find_gt_columns

        df = pd.DataFrame({"GENE": ["A"], "CHROM": ["1"]})
        assert _find_gt_columns(df) == []

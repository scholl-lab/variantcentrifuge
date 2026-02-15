"""Tests for GT column pre-parsing and caching.

Tests verify that GT column parsing happens once at DataFrame load time,
cache columns are properly cleaned up before output, and gene_burden.py
maintains no GT regex parsing (Phase 7 dead code removal).
"""

import re
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.dataframe_optimizer import (
    load_optimized_dataframe,
    parse_gt_column,
)


@pytest.mark.unit
class TestGTColumnParsing:
    """Tests for parse_gt_column function."""

    def test_parse_gt_column_basic(self):
        """Parse a simple GT string, verify output structure."""
        df = pd.DataFrame(
            {
                "GT": ["Sample1(0/1);Sample2(1/1)", "Sample3(0/0);Sample4(1/0)"],
                "CHROM": ["chr1", "chr2"],
            }
        )

        result_df = parse_gt_column(df)

        # Check column exists
        assert "_GT_PARSED" in result_df.columns

        # Check first row
        row0 = result_df["_GT_PARSED"][0]
        assert len(row0) == 2
        assert row0[0] == {"sample": "Sample1", "gt": "0/1"}
        assert row0[1] == {"sample": "Sample2", "gt": "1/1"}

        # Check second row
        row1 = result_df["_GT_PARSED"][1]
        assert len(row1) == 2
        assert row1[0] == {"sample": "Sample3", "gt": "0/0"}
        assert row1[1] == {"sample": "Sample4", "gt": "1/0"}

    def test_parse_gt_column_empty_values(self):
        """NaN, empty string, None all produce empty list."""
        df = pd.DataFrame(
            {
                "GT": [pd.NA, "", None, "   ", "Sample1(0/1)"],
                "CHROM": ["chr1", "chr2", "chr3", "chr4", "chr5"],
            }
        )

        result_df = parse_gt_column(df)

        # First four should be empty lists
        assert result_df["_GT_PARSED"][0] == []
        assert result_df["_GT_PARSED"][1] == []
        assert result_df["_GT_PARSED"][2] == []
        assert result_df["_GT_PARSED"][3] == []

        # Fifth should have one entry
        assert len(result_df["_GT_PARSED"][4]) == 1
        assert result_df["_GT_PARSED"][4][0] == {"sample": "Sample1", "gt": "0/1"}

    def test_parse_gt_column_single_sample(self):
        """Single sample GT string."""
        df = pd.DataFrame({"GT": ["Sample1(1/1)"], "CHROM": ["chr1"]})

        result_df = parse_gt_column(df)

        row0 = result_df["_GT_PARSED"][0]
        assert len(row0) == 1
        assert row0[0] == {"sample": "Sample1", "gt": "1/1"}

    def test_parse_gt_column_reference_genotype(self):
        """Includes 0/0 and ./. entries (still parsed, filtering is downstream)."""
        df = pd.DataFrame({"GT": ["Sample1(0/0);Sample2(./.)"], "CHROM": ["chr1"]})

        result_df = parse_gt_column(df)

        row0 = result_df["_GT_PARSED"][0]
        assert len(row0) == 2
        assert row0[0] == {"sample": "Sample1", "gt": "0/0"}
        assert row0[1] == {"sample": "Sample2", "gt": "./."}

    def test_parse_gt_column_no_gt_column(self):
        """DataFrame without GT column should skip parsing."""
        df = pd.DataFrame({"CHROM": ["chr1", "chr2"], "POS": [100, 200]})

        result_df = parse_gt_column(df)

        # Should not have _GT_PARSED column
        assert "_GT_PARSED" not in result_df.columns

    def test_parse_gt_column_empty_dataframe(self):
        """Empty DataFrame should skip parsing but add empty column."""
        df = pd.DataFrame({"GT": [], "CHROM": []})

        result_df = parse_gt_column(df)

        # Should have _GT_PARSED column but it should be empty
        assert "_GT_PARSED" in result_df.columns
        assert len(result_df["_GT_PARSED"]) == 0


@pytest.mark.unit
class TestGTCacheIntegration:
    """Tests for GT cache integration with load_optimized_dataframe."""

    def test_load_optimized_with_gt_cache(self, tmp_path):
        """Load TSV with GT column, verify _GT_PARSED column exists and is correct."""
        tsv_file = tmp_path / "test_variants.tsv"
        tsv_content = (
            "CHROM\tPOS\tREF\tALT\tGT\n"
            "chr1\t100\tA\tT\tSample1(0/1);Sample2(1/1)\n"
            "chr2\t200\tG\tC\tSample3(0/0)\n"
        )
        tsv_file.write_text(tsv_content)

        df, _rename_map = load_optimized_dataframe(str(tsv_file))

        # Check GT cache column exists
        assert "_GT_PARSED" in df.columns

        # Verify first row
        row0 = df["_GT_PARSED"][0]
        assert len(row0) == 2
        assert row0[0] == {"sample": "Sample1", "gt": "0/1"}
        assert row0[1] == {"sample": "Sample2", "gt": "1/1"}

        # Verify second row
        row1 = df["_GT_PARSED"][1]
        assert len(row1) == 1
        assert row1[0] == {"sample": "Sample3", "gt": "0/0"}

    def test_cache_columns_not_in_output(self):
        """Verify that columns starting with _ are dropped by cleanup logic."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2"],
                "POS": [100, 200],
                "_GT_PARSED": [
                    [{"sample": "S1", "gt": "0/1"}],
                    [{"sample": "S2", "gt": "1/1"}],
                ],
                "_INTERNAL_CACHE": ["val1", "val2"],
            }
        )

        # Simulate cache column cleanup (same logic as in output_stages.py)
        cache_cols = [c for c in df.columns if c.startswith("_")]
        assert len(cache_cols) == 2
        assert "_GT_PARSED" in cache_cols
        assert "_INTERNAL_CACHE" in cache_cols

        # Drop cache columns
        df_clean = df.drop(columns=cache_cols)

        # Verify cache columns are gone
        assert "_GT_PARSED" not in df_clean.columns
        assert "_INTERNAL_CACHE" not in df_clean.columns
        assert "CHROM" in df_clean.columns
        assert "POS" in df_clean.columns


@pytest.mark.unit
class TestGeneBurdenNoGTRegex:
    """Regression guard for Phase 7 dead GT loop removal in gene_burden.py."""

    def test_gene_burden_has_no_gt_regex(self):
        """Reads gene_burden.py source and asserts no re.compile calls with GT patterns exist.

        This is a regression guard ensuring Phase 7's dead code removal stays removed.
        gene_burden.py had a GT parsing loop that was never actually used (values were
        never read from the dict). Phase 7 removed it. This test ensures it doesn't
        creep back in.
        """
        # Find gene_burden.py
        base_path = Path(__file__).parent.parent.parent / "variantcentrifuge"
        gene_burden_path = base_path / "gene_burden.py"
        assert gene_burden_path.exists(), f"gene_burden.py not found at {gene_burden_path}"

        # Read source
        source = gene_burden_path.read_text()

        # Check for re.compile calls
        compile_pattern = re.compile(r"re\.compile\([^)]*\)")
        matches = compile_pattern.findall(source)

        # Should find zero re.compile calls (no regex in gene_burden.py)
        # If this fails, it means regex was added back - check if it's GT-related
        assert len(matches) == 0, (
            f"Found {len(matches)} re.compile call(s) in gene_burden.py: {matches}\n"
            "gene_burden.py should not have any regex compilation (dead GT loop removed in Phase 7)"
        )

        # Also check for explicit GT parsing patterns
        gt_patterns = [
            r"\([^()]+\)\([^)]+\)",  # The GT regex pattern
            "GT.*parse",  # Any GT parsing mention
            "genotype.*regex",  # Any genotype regex mention
        ]

        for pattern in gt_patterns:
            pattern_re = re.compile(pattern, re.IGNORECASE)
            matches = pattern_re.findall(source)
            assert len(matches) == 0, (
                f"Found GT-related pattern '{pattern}' in gene_burden.py: {matches}\n"
                "GT parsing should not exist in gene_burden.py (dead code removed in Phase 7)"
            )

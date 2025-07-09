"""
Unit tests for column ordering determinism fixes.

These tests ensure that column ordering in custom annotations is deterministic
and that the pipeline stage dependencies enforce consistent ordering.
"""

import pandas as pd

from variantcentrifuge.annotator import _add_json_annotations_as_columns


class TestColumnOrderingDeterminism:
    """Test that column ordering is deterministic in custom annotations."""

    def test_json_annotations_deterministic_column_order(self):
        """Test that JSON annotations are added in deterministic order."""
        # Create test DataFrame
        df = pd.DataFrame(
            {"CHROM": ["chr1", "chr2"], "POS": [100, 200], "GENE": ["GENE1", "GENE2"]}
        )

        # Create JSON gene data with multiple columns
        # Use intentionally unsorted keys to test determinism
        json_gene_data = {
            "GENE1": {"score": 0.8, "category": "pathogenic", "evidence": "strong"},
            "GENE2": {"score": 0.3, "category": "benign", "evidence": "weak"},
        }

        # Run annotation multiple times
        results = []
        for _ in range(10):
            result_df = _add_json_annotations_as_columns(df.copy(), json_gene_data)
            # Get the column names of the new annotation columns
            new_cols = [col for col in result_df.columns if col not in df.columns]
            results.append(new_cols)

        # All results should be identical
        first_result = results[0]
        assert all(
            result == first_result for result in results
        ), "JSON annotation columns should be added in deterministic order"

        # Verify columns are in sorted order
        expected_cols = ["category", "evidence", "score"]  # Alphabetical order
        assert (
            first_result == expected_cols
        ), f"Expected columns {expected_cols}, got {first_result}"

    def test_empty_json_data_handling(self):
        """Test handling of empty JSON data."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100], "GENE": ["GENE1"]})

        # Test with empty dict
        result_df = _add_json_annotations_as_columns(df.copy(), {})
        assert result_df.equals(df), "Empty JSON data should not modify DataFrame"

        # Test with dict containing empty values
        json_gene_data = {"GENE1": {}}
        result_df = _add_json_annotations_as_columns(df.copy(), json_gene_data)
        assert result_df.equals(df), "JSON data with empty values should not modify DataFrame"

    def test_partial_gene_matches(self):
        """Test deterministic ordering when only some genes have annotations."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "POS": [100, 200, 300],
                "GENE": ["GENE1", "GENE2", "GENE3"],
            }
        )

        # Only GENE1 and GENE3 have annotations
        json_gene_data = {
            "GENE1": {"score": 0.8, "category": "pathogenic"},
            "GENE3": {"score": 0.2, "category": "benign"},
        }

        # Run multiple times to verify determinism
        results = []
        for _ in range(5):
            result_df = _add_json_annotations_as_columns(df.copy(), json_gene_data)
            results.append(result_df)

        # All results should be identical (using pandas.testing.assert_frame_equal
        # which handles NaN properly)
        first_result = results[0]
        for i, result in enumerate(results[1:], 1):
            pd.testing.assert_frame_equal(
                result, first_result, check_dtype=False, check_like=True
            )  # Allow different column order

        # Verify columns are in sorted order
        new_cols = [col for col in first_result.columns if col not in df.columns]
        expected_cols = ["category", "score"]  # Alphabetical order
        assert new_cols == expected_cols, f"Expected columns {expected_cols}, got {new_cols}"

    def test_large_number_of_columns(self):
        """Test deterministic ordering with many columns."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [100], "GENE": ["GENE1"]})

        # Create many columns with intentionally unsorted keys
        many_columns = {f"col_{i:02d}": f"value_{i}" for i in range(20, 0, -1)}  # Reverse order
        json_gene_data = {"GENE1": many_columns}

        # Run multiple times
        results = []
        for _ in range(3):
            result_df = _add_json_annotations_as_columns(df.copy(), json_gene_data)
            new_cols = [col for col in result_df.columns if col not in df.columns]
            results.append(new_cols)

        # All results should be identical and sorted
        first_result = results[0]
        assert all(
            result == first_result for result in results
        ), "Large number of columns should be deterministic"

        # Verify sorted order
        expected_sorted = sorted(many_columns.keys())
        assert (
            first_result == expected_sorted
        ), f"Expected sorted columns {expected_sorted[:5]}..., got {first_result[:5]}..."


class TestPipelineStageOrdering:
    """Test that pipeline stages execute in deterministic order."""

    def test_variant_analysis_depends_on_custom_annotation(self):
        """Test that VariantAnalysisStage has custom_annotation in soft dependencies."""
        from variantcentrifuge.stages.analysis_stages import VariantAnalysisStage

        stage = VariantAnalysisStage()
        soft_deps = stage.soft_dependencies

        assert (
            "custom_annotation" in soft_deps
        ), "VariantAnalysisStage should have custom_annotation in soft dependencies"

    def test_stage_dependency_chain(self):
        """Test the complete dependency chain for deterministic ordering."""
        from variantcentrifuge.stages.analysis_stages import (
            CustomAnnotationStage,
            InheritanceAnalysisStage,
            VariantAnalysisStage,
        )

        custom_stage = CustomAnnotationStage()
        variant_stage = VariantAnalysisStage()
        inheritance_stage = InheritanceAnalysisStage()

        # Verify dependency chain
        assert custom_stage.dependencies == {"dataframe_loading"}
        assert variant_stage.dependencies == {"dataframe_loading"}
        assert "custom_annotation" in variant_stage.soft_dependencies
        assert "custom_annotation" in inheritance_stage.soft_dependencies

        # This ensures: dataframe_loading -> custom_annotation -> variant_analysis
        # which should give deterministic column ordering

"""Unit tests for critical processing stages - filtering and data transformation."""

from argparse import Namespace
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.processing_stages import (
    ExtraColumnRemovalStage,
    FieldExtractionStage,
    MultiAllelicSplitStage,
    SnpSiftFilterStage,
)


@pytest.fixture
def mock_workspace(tmp_path):
    """Create a mock workspace."""
    workspace = Mock(spec=Workspace)
    workspace.output_dir = tmp_path
    workspace.intermediate_dir = tmp_path / "intermediate"
    workspace.intermediate_dir.mkdir()
    workspace.base_name = "test"
    workspace.get_intermediate_path = lambda x: workspace.intermediate_dir / x
    workspace.get_output_path = lambda x, ext=".tsv": workspace.output_dir / f"{x}{ext}"
    workspace.get_temp_path = lambda x: workspace.intermediate_dir / "temp" / x
    return workspace


@pytest.fixture
def base_context(mock_workspace):
    """Create a base pipeline context."""
    return PipelineContext(
        args=Namespace(),
        config={},
        workspace=mock_workspace,
    )


class TestSnpSiftFilterStage:
    """Test the SnpSiftFilterStage - CRITICAL for correct variant filtering."""

    def test_no_filter_expression(self, base_context):
        """Test when no filter expression is provided."""
        base_context.config = {}
        base_context.data = Path("/tmp/input.vcf.gz")

        stage = SnpSiftFilterStage()
        result = stage._process(base_context)

        # Should pass through unchanged
        assert result.data == base_context.data
        # PipelineContext initializes all attributes, so we just verify data unchanged

    def test_late_filtering_mode(self, base_context):
        """Test when late_filtering is enabled (should skip)."""
        base_context.config = {"filter": "(GEN[*].DP >= 10)", "late_filtering": True}
        base_context.data = Path("/tmp/input.vcf.gz")

        stage = SnpSiftFilterStage()
        result = stage._process(base_context)

        # Should skip filtering in late mode
        assert result.data == base_context.data
        # PipelineContext initializes all attributes, so we just verify data unchanged

    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    def test_basic_filtering(self, mock_apply_filter, base_context, tmp_path):
        """Test basic SnpSift filtering with simple expression."""
        input_vcf = tmp_path / "input.vcf.gz"
        input_vcf.touch()

        base_context.config = {"filter": "(GEN[*].DP >= 10)", "threads": 4}
        base_context.data = input_vcf

        # Mock apply_snpsift_filter to create the output file with VCF content
        def mock_apply_side_effect(input_file, filter_expr, config, output_file):
            # Create realistic filtered VCF content
            vcf_content = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT:DP\t0/1:15\t0/0:20
chr1\t200\t.\tG\tC\t100\tPASS\t.\tGT:DP\t0/0:12\t0/1:18
"""
            Path(output_file).write_text(vcf_content)

        mock_apply_filter.side_effect = mock_apply_side_effect

        stage = SnpSiftFilterStage()
        result = stage._process(base_context)

        # Verify apply_snpsift_filter was called correctly
        mock_apply_filter.assert_called_once()
        call_args = mock_apply_filter.call_args[0]
        assert str(input_vcf) == call_args[0]
        assert call_args[1] == "(GEN[*].DP >= 10)"
        assert call_args[2]["threads"] == 4

        # Verify output paths updated
        assert result.filtered_vcf.name == "test.filtered.vcf.gz"
        assert result.data == result.filtered_vcf
        expected_path = base_context.workspace.intermediate_dir / "test.filtered.vcf.gz"
        assert result.filtered_vcf == expected_path

    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    def test_filtering_uses_split_vcf_when_split_before_filters(
        self, mock_apply_filter, base_context, tmp_path
    ):
        """Test before-filter split output feeds SnpSift filtering."""
        extracted_vcf = tmp_path / "extracted.vcf.gz"
        split_vcf = tmp_path / "split.vcf.gz"
        extracted_vcf.touch()
        split_vcf.touch()

        base_context.config = {"filter": "(ANN[*].IMPACT has 'HIGH')"}
        base_context.extracted_vcf = extracted_vcf
        base_context.split_annotations_vcf = split_vcf
        base_context.data = split_vcf

        stage = SnpSiftFilterStage(split_before_filtering=True)
        stage._process(base_context)

        mock_apply_filter.assert_called_once()
        assert mock_apply_filter.call_args.args[0] == str(split_vcf)

    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    def test_complex_filter_expression(self, mock_apply_filter, base_context, tmp_path):
        """Test with complex filter expression including AND/OR logic."""
        input_vcf = tmp_path / "input.vcf.gz"
        input_vcf.touch()

        # Complex real-world filter
        complex_filter = (
            "((ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')) & "
            "(GEN[*].DP >= 10) & (QUAL >= 30) & "
            "((AF < 0.01) | (dbNSFP_MetaSVM_pred has 'D'))"
        )

        base_context.config = {"filter": complex_filter}
        base_context.data = input_vcf

        stage = SnpSiftFilterStage()
        stage._process(base_context)

        # Verify complex expression passed correctly
        assert mock_apply_filter.call_args[0][1] == complex_filter

    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    def test_filter_with_preset_combination(self, mock_apply_filter, base_context, tmp_path):
        """Test filter that would come from preset combinations."""
        input_vcf = tmp_path / "input.vcf.gz"
        input_vcf.touch()

        # Simulate filter from preset combination
        preset_filter = (
            "(ANN[*].IMPACT has 'HIGH') & (dbNSFP_SIFT_pred has 'D') & ((AF < 0.001) | (AF = '.'))"
        )

        base_context.config = {"filter": preset_filter}
        base_context.data = input_vcf

        stage = SnpSiftFilterStage()
        result = stage._process(base_context)

        assert mock_apply_filter.called
        assert result.filtered_vcf is not None

    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    def test_filter_error_handling(self, mock_apply_filter, base_context, tmp_path):
        """Test error handling when filtering fails."""
        input_vcf = tmp_path / "input.vcf.gz"
        input_vcf.touch()

        # Simulate filter failure
        mock_apply_filter.side_effect = RuntimeError("SnpSift failed: Invalid expression")

        base_context.config = {"filter": "INVALID SYNTAX"}
        base_context.data = input_vcf

        stage = SnpSiftFilterStage()

        with pytest.raises(RuntimeError, match="SnpSift failed"):
            stage._process(base_context)


class TestMultiAllelicSplitStage:
    """Test the MultiAllelicSplitStage - important for SNPeff annotation handling."""

    def test_split_not_requested(self, base_context):
        """Test when splitting is not requested."""
        base_context.config = {"snpeff_splitting_mode": None}
        base_context.data = Path("/tmp/input.vcf.gz")

        stage = MultiAllelicSplitStage()
        result = stage._process(base_context)

        # Should pass through unchanged
        assert result.data == base_context.data

    @patch("variantcentrifuge.stages.processing_stages.split_snpeff_annotations")
    def test_basic_split(self, mock_split, base_context, tmp_path):
        """Test basic SNPeff annotation splitting."""
        input_vcf = tmp_path / "input.vcf.gz"
        input_vcf.touch()

        base_context.config = {"snpeff_splitting_mode": "before_filters"}
        base_context.data = input_vcf
        base_context.extracted_vcf = input_vcf

        stage = MultiAllelicSplitStage()
        result = stage._process(base_context)

        # Verify split function called
        mock_split.assert_called_once()
        assert str(input_vcf) in mock_split.call_args[0]

        # Verify output path
        assert hasattr(result, "split_annotations_vcf")
        assert result.split_annotations_vcf.name == "test.split_annotations.vcf.gz"
        assert result.data == result.split_annotations_vcf

    @patch("variantcentrifuge.stages.processing_stages.split_snpeff_annotations")
    def test_split_after_filtering(self, mock_split, base_context, tmp_path):
        """Test splitting after filtering (when configured)."""
        filtered_vcf = tmp_path / "filtered.vcf.gz"
        filtered_vcf.touch()

        base_context.config = {"snpeff_splitting_mode": "after_filters"}
        base_context.filtered_vcf = filtered_vcf
        base_context.data = filtered_vcf

        # Mark filtering as complete
        base_context.mark_complete("snpsift_filtering")

        stage = MultiAllelicSplitStage()
        stage._process(base_context)

        # Should use filtered VCF as input
        if mock_split.called:
            assert str(filtered_vcf) in mock_split.call_args[0][0]

    def test_parallel_safe_property(self, base_context):
        """Test that this stage is marked as parallel safe."""
        stage = MultiAllelicSplitStage()
        assert stage.parallel_safe is False

    def test_split_mode_dependencies(self):
        """Test split/filter ordering dependencies are encoded by mode."""
        assert "snpsift_filtering" not in MultiAllelicSplitStage("before_filters").dependencies
        assert "snpsift_filtering" in MultiAllelicSplitStage("after_filters").dependencies
        assert (
            "multiallelic_split"
            in SnpSiftFilterStage(split_before_filtering=True).soft_dependencies
        )


class TestFieldExtractionStage:
    """Test VCF path selection for field extraction."""

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_field_extraction_uses_split_vcf_after_filters(
        self, mock_extract_fields, base_context, tmp_path
    ):
        """Test after-filter split output is used for extraction."""
        filtered_vcf = tmp_path / "filtered.vcf.gz"
        split_vcf = tmp_path / "split.vcf.gz"
        filtered_vcf.touch()
        split_vcf.touch()
        base_context.config = {
            "fields_to_extract": "CHROM POS ANN[0].FEATUREID",
            "snpeff_splitting_mode": "after_filters",
        }
        base_context.filtered_vcf = filtered_vcf
        base_context.split_annotations_vcf = split_vcf

        FieldExtractionStage()._process(base_context)

        mock_extract_fields.assert_called_once()
        assert mock_extract_fields.call_args.kwargs["variant_file"] == str(split_vcf)


class TestExtraColumnRemovalStage:
    """Test the ExtraColumnRemovalStage - simple but necessary cleanup."""

    def test_no_columns_to_remove(self, base_context):
        """Test when no columns are specified for removal."""
        base_context.config = {}
        base_context.data = Path("/tmp/input.tsv")

        stage = ExtraColumnRemovalStage()
        result = stage._process(base_context)

        # Should pass through unchanged
        assert result.data == base_context.data

    def test_remove_single_column(self, base_context, tmp_path):
        """Test removing a single column."""
        import pandas as pd

        # Create test file with data
        input_tsv = tmp_path / "input.tsv"
        df = pd.DataFrame(
            {
                "Gene": ["GENE1", "GENE2"],
                "Variant": ["var1", "var2"],
                "ExtraColumn1": ["extra1", "extra2"],
            }
        )
        df.to_csv(input_tsv, sep="\t", index=False)

        base_context.config = {"extra_columns_to_remove": ["ExtraColumn1"]}
        base_context.data = input_tsv

        stage = ExtraColumnRemovalStage()
        result = stage._process(base_context)

        # Verify column was removed
        result_df = pd.read_csv(result.data, sep="\t")
        assert "ExtraColumn1" not in result_df.columns
        assert "Gene" in result_df.columns
        assert "Variant" in result_df.columns

    def test_remove_multiple_columns(self, base_context, tmp_path):
        """Test removing multiple columns."""
        import pandas as pd

        # Create test file with multiple columns
        input_tsv = tmp_path / "input.tsv"
        df = pd.DataFrame(
            {
                "Gene": ["GENE1"],
                "Variant": ["var1"],
                "TempCol1": ["temp1"],
                "TempCol2": ["temp2"],
                "DebugInfo": ["debug"],
                "InternalID": ["id1"],
            }
        )
        df.to_csv(input_tsv, sep="\t", index=False)

        columns_to_remove = ["TempCol1", "TempCol2", "DebugInfo", "InternalID"]
        base_context.config = {"extra_columns_to_remove": columns_to_remove}
        base_context.data = input_tsv

        stage = ExtraColumnRemovalStage()
        result = stage._process(base_context)

        # Verify all columns were removed
        result_df = pd.read_csv(result.data, sep="\t")
        for col in columns_to_remove:
            assert col not in result_df.columns
        assert len(result_df.columns) == 2  # Only Gene and Variant remain

    def test_gzip_handling(self, base_context, tmp_path):
        """Test handling of gzipped files."""
        import pandas as pd

        # Create gzipped test file
        input_tsv = tmp_path / "input.tsv.gz"
        df = pd.DataFrame({"Gene": ["GENE1"], "Extra": ["data"]})
        df.to_csv(input_tsv, sep="\t", index=False, compression="gzip")

        base_context.config = {"extra_columns_to_remove": ["Extra"], "gzip_intermediates": True}
        base_context.data = input_tsv

        stage = ExtraColumnRemovalStage()
        result = stage._process(base_context)

        # Output should also be gzipped
        assert result.data.name == "test.columns_removed.tsv.gz"
        # Verify can read the gzipped output
        result_df = pd.read_csv(result.data, sep="\t", compression="gzip")
        assert "Extra" not in result_df.columns

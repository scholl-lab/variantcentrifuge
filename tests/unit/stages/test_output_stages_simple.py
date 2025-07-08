"""Unit tests for output stages - simplified version."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch
from concurrent.futures import Future
import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    VariantIdentifierStage,
    FinalFilteringStage,
    TSVOutputStage,
    MetadataGenerationStage,
    ArchiveCreationStage,
    ParallelReportGenerationStage,
    ExcelReportStage,
    HTMLReportStage,
)
from tests.mocks.fixtures import create_test_context


class TestVariantIdentifierStage:
    """Test VariantIdentifierStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        # Create a test DataFrame
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2"],
                "POS": [100, 200, 300],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "G"],
                "QUAL": [30, 40, 50],
            }
        )
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_add_variant_identifier(self, context):
        """Test adding variant identifiers."""
        stage = VariantIdentifierStage()
        result = stage(context)

        # Check that VAR_ID column was added (implementation uses VAR_ID, not Variant_ID)
        assert "VAR_ID" in result.current_dataframe.columns

        # Check format of variant IDs - implementation uses var_XXXX_hash format
        variant_ids = result.current_dataframe["VAR_ID"].tolist()
        assert variant_ids[0].startswith("var_0001_")
        assert variant_ids[1].startswith("var_0002_")
        assert variant_ids[2].startswith("var_0003_")

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = VariantIdentifierStage()
        assert stage.dependencies == {"dataframe_loading"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = VariantIdentifierStage()
        # Implementation marks this as False because it modifies DataFrame in-place
        assert stage.parallel_safe is False


class TestFinalFilteringStage:
    """Test FinalFilteringStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        ctx.current_dataframe = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr1", "chr2", "chr3"],
                "POS": [100, 200, 300, 400],
                "QUAL": [30, 40, 20, 50],
                "AF": [0.01, 0.1, 0.001, 0.5],
            }
        )
        ctx.mark_complete("dataframe_loading")
        return ctx

    def test_final_filtering(self, context):
        """Test final filtering with pandas query."""
        context.config["final_filter"] = "QUAL >= 30 & AF < 0.1"

        stage = FinalFilteringStage()
        result = stage(context)

        # Check filtering worked (only chr1:100 should remain)
        assert len(result.current_dataframe) == 1
        assert all(result.current_dataframe["QUAL"] >= 30)
        assert all(result.current_dataframe["AF"] < 0.1)

    def test_no_filtering(self, context):
        """Test when no filtering is needed."""
        # No filters configured
        stage = FinalFilteringStage()
        result = stage(context)

        # DataFrame should be unchanged
        assert len(result.current_dataframe) == 4

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = FinalFilteringStage()
        assert stage.parallel_safe is True


class TestTSVOutputStage:
    """Test TSVOutputStage."""

    @pytest.fixture
    def context(self):
        """Create test context with DataFrame."""
        ctx = create_test_context()
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.current_dataframe = pd.DataFrame(
            {"CHROM": ["chr1", "chr1"], "POS": [100, 200], "QUAL": [30, 40]}
        )
        ctx.mark_complete("dataframe_loading")
        return ctx

    @patch("pandas.DataFrame.to_csv")
    def test_tsv_output(self, mock_to_csv, context):
        """Test TSV file output."""
        # TSVOutputStage depends on variant_identifier, so mark it complete
        context.mark_complete("variant_identifier")

        stage = TSVOutputStage()
        result = stage(context)

        # Check to_csv was called
        assert mock_to_csv.called

        # Check context updated
        assert hasattr(result, "final_output_file") or hasattr(result, "final_output_path")

    def test_missing_dataframe_error(self, context):
        """Test error when no DataFrame available."""
        context.current_dataframe = None
        # TSVOutputStage depends on variant_identifier, so mark it complete
        context.mark_complete("variant_identifier")

        stage = TSVOutputStage()
        # Should log warning and return
        result = stage(context)
        assert result == context


class TestMetadataGenerationStage:
    """Test MetadataGenerationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(gene_name="BRCA1", vcf_file="/tmp/input.vcf")
        ctx.config["output_file"] = "/tmp/output.tsv"
        ctx.mark_complete("tsv_output")
        return ctx

    @patch("variantcentrifuge.stages.output_stages.sanitize_metadata_field")
    @patch("variantcentrifuge.stages.output_stages.get_tool_version")
    def test_metadata_generation(self, mock_get_version, mock_sanitize, context):
        """Test metadata file generation."""
        # Add normalized_genes to config
        context.config["normalized_genes"] = ["BRCA1"]

        # Mock tool versions
        mock_get_version.return_value = "1.0.0"

        # Mock sanitize to pass through the metadata unchanged
        mock_sanitize.side_effect = lambda x: x

        stage = MetadataGenerationStage()
        result = stage(context)

        # Check metadata file was created
        assert "metadata" in result.report_paths
        metadata_path = result.report_paths["metadata"]
        assert metadata_path.exists()
        assert metadata_path.suffix == ".json"

        # Read and check contents
        with open(metadata_path) as f:
            metadata = json.load(f)

        # Check metadata contents
        assert metadata["input_files"]["genes"] == ["BRCA1"]
        assert metadata["input_files"]["vcf"] == "/tmp/input.vcf"
        assert "run_date" in metadata

    def test_stage_properties(self):
        """Test stage properties."""
        stage = MetadataGenerationStage()
        assert stage.dependencies == {"tsv_output"}
        assert stage.parallel_safe is True


class TestArchiveCreationStage:
    """Test ArchiveCreationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context(config_overrides={"archive_results": True})
        # Create a real temp directory
        ctx.workspace.output_dir = Path(tempfile.mkdtemp())
        ctx.mark_complete("tsv_output")
        return ctx

    def test_archive_creation(self, context):
        """Test archive creation."""
        # Create a simple file to archive
        test_file = context.workspace.output_dir / "test.txt"
        test_file.write_text("test content")

        stage = ArchiveCreationStage()
        result = stage(context)

        # Check archive was created
        assert "archive" in result.report_paths
        archive_path = result.report_paths["archive"]
        assert archive_path.exists()
        assert archive_path.suffix == ".gz"

    def test_skip_if_disabled(self, context):
        """Test skipping when archiving disabled."""
        context.config["archive_results"] = False

        stage = ArchiveCreationStage()
        result = stage(context)

        # Should return unchanged context
        assert result == context

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = ArchiveCreationStage()
        assert stage.parallel_safe is True


class TestParallelReportGenerationStage:
    """Test ParallelReportGenerationStage."""

    @pytest.fixture
    def context(self):
        """Create test context."""
        ctx = create_test_context()
        ctx.data = Path("/tmp/test.tsv")
        ctx.mark_complete("tsv_output")
        return ctx

    def test_no_reports_requested(self, context):
        """Test when no reports are requested."""
        context.config = {"no_metadata": True}

        stage = ParallelReportGenerationStage()
        result = stage(context)

        # Should pass through unchanged
        assert result == context

    @patch.object(MetadataGenerationStage, "_process")
    def test_single_report_generation(self, mock_metadata, context):
        """Test generating a single report."""
        context.config = {"no_metadata": False}
        mock_metadata.return_value = context

        stage = ParallelReportGenerationStage()
        stage(context)

        # Should call the single report stage
        assert mock_metadata.called

    @patch("variantcentrifuge.stages.output_stages.as_completed")
    @patch("variantcentrifuge.stages.output_stages.ThreadPoolExecutor")
    @patch.object(ExcelReportStage, "_process")
    @patch.object(HTMLReportStage, "_process")
    @patch.object(MetadataGenerationStage, "_process")
    def test_parallel_report_generation(
        self, mock_metadata, mock_html, mock_excel, mock_executor, mock_as_completed, context
    ):
        """Test generating multiple reports in parallel."""
        context.config = {"xlsx": True, "html_report": True, "no_metadata": False}

        # Mock the report stages
        mock_excel.return_value = context
        mock_html.return_value = context
        mock_metadata.return_value = context

        # Mock futures
        mock_futures = []
        for _ in range(3):
            future = Mock(spec=Future)
            future.result.return_value = None
            mock_futures.append(future)

        # Mock the executor
        mock_executor_instance = Mock()
        mock_executor_instance.submit.side_effect = mock_futures
        mock_executor.return_value.__enter__.return_value = mock_executor_instance

        # Mock as_completed to return futures immediately
        mock_as_completed.return_value = mock_futures

        stage = ParallelReportGenerationStage()
        stage(context)

        # Should use ThreadPoolExecutor
        assert mock_executor.called
        # Should submit 3 tasks (excel, html, metadata)
        assert mock_executor_instance.submit.call_count == 3

    @patch("variantcentrifuge.stages.output_stages.as_completed")
    @patch("variantcentrifuge.stages.output_stages.ThreadPoolExecutor")
    @patch.object(ExcelReportStage, "_process")
    def test_report_generation_error_handling(
        self, mock_excel, mock_executor, mock_as_completed, context
    ):
        """Test error handling during report generation."""
        context.config = {"excel": True}

        # Mock excel generation to fail
        mock_excel.side_effect = Exception("Excel generation failed")

        # Mock the executor
        mock_future = Mock(spec=Future)
        mock_future.result.side_effect = Exception("Excel generation failed")
        mock_executor_instance = Mock()
        mock_executor_instance.submit.return_value = mock_future
        mock_executor.return_value.__enter__.return_value = mock_executor_instance

        # Mock as_completed
        mock_as_completed.return_value = [mock_future]

        stage = ParallelReportGenerationStage()
        # Should not raise exception, just log error
        result = stage(context)

        assert result == context

    def test_parallel_safe_property(self, context):
        """Test that this stage is not marked as parallel safe."""
        stage = ParallelReportGenerationStage()
        assert stage.parallel_safe is False  # Manages its own parallelism

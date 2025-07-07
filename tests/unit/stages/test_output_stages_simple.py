"""Unit tests for output stages - simplified version."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd
import pytest

from variantcentrifuge.stages.output_stages import (
    VariantIdentifierStage,
    FinalFilteringStage,
    TSVOutputStage,
    MetadataGenerationStage,
    ArchiveCreationStage,
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

        # Check that Variant_ID column was added
        assert "Variant_ID" in result.current_dataframe.columns

        # Check format of variant IDs
        variant_ids = result.current_dataframe["Variant_ID"].tolist()
        assert variant_ids[0] == "chr1:100:A>T"
        assert variant_ids[1] == "chr1:200:G>C"
        assert variant_ids[2] == "chr2:300:C>G"

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = VariantIdentifierStage()
        assert stage.dependencies == {"dataframe_loading"}

    def test_parallel_safe(self):
        """Test parallel safety."""
        stage = VariantIdentifierStage()
        assert stage.parallel_safe is True


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
        stage = TSVOutputStage()
        result = stage(context)

        # Check to_csv was called
        assert mock_to_csv.called

        # Check context updated
        assert hasattr(result, "final_output_file") or hasattr(result, "final_output_path")

    def test_missing_dataframe_error(self, context):
        """Test error when no DataFrame available."""
        context.current_dataframe = None

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

    def test_metadata_generation(self, context):
        """Test metadata file generation."""
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
        assert metadata["gene_name"] == "BRCA1"
        assert metadata["vcf_file"] == "/tmp/input.vcf"
        assert "timestamp" in metadata

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

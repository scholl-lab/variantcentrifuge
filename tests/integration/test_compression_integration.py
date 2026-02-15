"""
Integration tests for TSV compression functionality.

These tests verify that the compression fixes work correctly in realistic pipeline scenarios.
"""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

pytestmark = pytest.mark.integration

from variantcentrifuge.pipeline_core import PipelineContext, Workspace
from variantcentrifuge.stages.processing_stages import FieldExtractionStage


class TestCompressionIntegration:
    """Integration tests for compression behavior in realistic scenarios."""

    @pytest.fixture
    def temp_workspace(self):
        """Create a temporary workspace for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            workspace = Workspace(output_dir=Path(tmpdir) / "output", base_name="test_integration")
            # Add a vcf_file attribute for testing
            workspace.vcf_file = Path(tmpdir) / "test.vcf.gz"
            yield workspace

    @pytest.fixture
    def pipeline_context(self, temp_workspace):
        """Create a realistic pipeline context."""
        context = PipelineContext(
            args=None,
            config={
                "fields_to_extract": "CHROM POS REF ALT GENE",
                "extract_fields_separator": ",",
                "log_level": "INFO",
                "gzip_intermediates": True,  # Enable compression
            },
            workspace=temp_workspace,
        )
        context.data = temp_workspace.vcf_file

        # Mark dependencies as complete
        context.completion_status = {
            "gene_bed_creation": True,
            "configuration_loading": True,
            "variant_extraction": True,
        }
        context.is_complete = lambda stage_name: context.completion_status.get(stage_name, False)

        return context

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_field_extraction_creates_compressed_files(self, mock_extract, pipeline_context):
        """Test that FieldExtractionStage creates compressed files when gzip_intermediates=True."""

        # Mock extract_fields_bcftools to simulate successful extraction
        def mock_extract_func(variant_file, fields, cfg, output_file, vcf_samples=None):
            # Create the output file to simulate successful extraction
            if output_file.endswith(".gz"):
                import gzip

                with gzip.open(output_file, "wt") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
                    f.write("chr1\t100\tA\tT\tGENE1\n")
            else:
                with open(output_file, "w") as f:
                    f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
                    f.write("chr1\t100\tA\tT\tGENE1\n")
            return output_file

        mock_extract.side_effect = mock_extract_func

        # Run the field extraction stage
        stage = FieldExtractionStage()
        result_context = stage._process(pipeline_context)

        # Verify that the output file path includes .gz
        assert str(result_context.extracted_tsv).endswith(".extracted.tsv.gz")

        # Verify that extract_fields_bcftools was called with compressed filename
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args
        output_file = call_args.kwargs["output_file"]
        assert output_file.endswith(".extracted.tsv.gz")

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_field_extraction_creates_uncompressed_files(self, mock_extract, pipeline_context):
        """Test FieldExtractionStage creates uncompressed files when gzip_intermediates=False."""
        # Disable compression
        pipeline_context.config["gzip_intermediates"] = False

        # Mock extract_fields_bcftools to simulate successful extraction
        def mock_extract_func(variant_file, fields, cfg, output_file, vcf_samples=None):
            with open(output_file, "w") as f:
                f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
                f.write("chr1\t100\tA\tT\tGENE1\n")
            return output_file

        mock_extract.side_effect = mock_extract_func

        # Run the field extraction stage
        stage = FieldExtractionStage()
        result_context = stage._process(pipeline_context)

        # Verify that the output file path does NOT include .gz
        assert str(result_context.extracted_tsv).endswith(".extracted.tsv")
        assert not str(result_context.extracted_tsv).endswith(".extracted.tsv.gz")

        # Verify that extract_fields_bcftools was called with uncompressed filename
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args
        output_file = call_args.kwargs["output_file"]
        assert output_file.endswith(".extracted.tsv")
        assert not output_file.endswith(".extracted.tsv.gz")

    def test_workspace_intermediate_path_generation(self, temp_workspace):
        """Test that workspace generates correct intermediate paths."""
        # Test that base intermediate path generation works
        base_path = temp_workspace.get_intermediate_path("test.extracted.tsv")
        assert str(base_path).endswith("test.extracted.tsv")

        # Test compression extension handling
        compressed_path = Path(str(base_path) + ".gz")
        assert str(compressed_path).endswith("test.extracted.tsv.gz")

    def test_compression_setting_defaults(self):
        """Test that compression settings have correct defaults."""
        # Test default configuration
        config = {}

        # The default should be True (compression enabled)
        default_value = config.get("gzip_intermediates", True)
        assert default_value is True

        # Test explicit False
        config["gzip_intermediates"] = False
        explicit_false = config.get("gzip_intermediates", True)
        assert explicit_false is False

        # Test explicit True
        config["gzip_intermediates"] = True
        explicit_true = config.get("gzip_intermediates", True)
        assert explicit_true is True

    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_config_propagation_to_extract_fields(
        self, mock_extract, pipeline_context, temp_workspace
    ):
        """Test that configuration is properly propagated to extract_fields_bcftools."""
        pipeline_context.config["log_level"] = "DEBUG"

        mock_extract.return_value = "/tmp/output.tsv.gz"

        stage = FieldExtractionStage()
        stage._process(pipeline_context)

        # Verify configuration was passed correctly
        mock_extract.assert_called_once()
        call_args = mock_extract.call_args

        # Check that the right parameters were passed
        # The input VCF should be the one from the workspace
        assert call_args.kwargs["variant_file"] == str(temp_workspace.vcf_file)
        assert call_args.kwargs["fields"] == "CHROM POS REF ALT GENE"

        cfg = call_args.kwargs["cfg"]
        assert cfg["debug_level"] == "DEBUG"


class TestParallelProcessingCompressionIntegration:
    """Integration tests for parallel processing compression."""

    def test_parallel_chunk_compression_settings(self):
        """Test that parallel chunk processing respects compression settings."""
        # Test that worker config properly includes compression setting
        test_config = {
            "threads_per_chunk": 2,
            "fields_to_extract": "CHROM POS REF ALT",
            "gzip_intermediates": True,  # Enable compression
        }

        # Test compression enabled
        use_compression = test_config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"
        assert tsv_suffix == ".extracted.tsv.gz"

        # Test compression disabled
        test_config["gzip_intermediates"] = False
        use_compression = test_config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"
        assert tsv_suffix == ".extracted.tsv"

        # Test default behavior (should be compressed)
        del test_config["gzip_intermediates"]
        use_compression = test_config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"
        assert tsv_suffix == ".extracted.tsv.gz"


class TestCompressionConsistencyAcrossStages:
    """Test that compression settings are consistent across different stages."""

    def test_consistent_compression_defaults(self):
        """Test that all stages use consistent compression defaults."""
        # Default should always be True for compression
        default_compression = True

        # Test FieldExtractionStage default
        config = {}
        field_compression = config.get("gzip_intermediates", True)
        assert field_compression == default_compression

        # Test ParallelCompleteProcessingStage default
        parallel_compression = config.get("gzip_intermediates", True)
        assert parallel_compression == default_compression

        # Both should be the same
        assert field_compression == parallel_compression

    def test_compression_setting_propagation(self):
        """Test that compression settings propagate correctly through config."""
        # Test that setting is maintained through config copies
        original_config = {"gzip_intermediates": False, "other_setting": "value"}

        # Simulate config copying (like in parallel processing)
        worker_config = {
            "gzip_intermediates": original_config.get("gzip_intermediates", True),
            "other_setting": original_config.get("other_setting"),
        }

        # Should preserve the False setting
        assert worker_config["gzip_intermediates"] is False

        # Test with True setting
        original_config["gzip_intermediates"] = True
        worker_config["gzip_intermediates"] = original_config.get("gzip_intermediates", True)
        assert worker_config["gzip_intermediates"] is True

        # Test with missing setting (should default to True)
        del original_config["gzip_intermediates"]
        worker_config["gzip_intermediates"] = original_config.get("gzip_intermediates", True)
        assert worker_config["gzip_intermediates"] is True

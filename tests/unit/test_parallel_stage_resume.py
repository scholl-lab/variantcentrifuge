"""Tests for parallel stage resume functionality."""

import gzip
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.stages.processing_stages import (
    ParallelCompleteProcessingStage,
    ParallelVariantExtractionStage,
)


class TestParallelVariantExtractionStage:
    """Test parallel variant extraction stage resume functionality."""

    def test_validate_existing_chunk_valid_vcf(self):
        """Test validation of valid VCF chunk."""
        stage = ParallelVariantExtractionStage()

        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tmp:
            # Create a valid VCF file
            with gzip.open(tmp.name, "wt") as f:
                f.write("##fileformat=VCFv4.2\n")
                f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                f.write("chr1\t1000\t.\tA\tT\t30\tPASS\t.\n")

            chunk_path = Path(tmp.name)
            assert stage._validate_existing_chunk(chunk_path) is True

            # Cleanup
            chunk_path.unlink()

    def test_validate_existing_chunk_empty_file(self):
        """Test validation rejects empty files."""
        stage = ParallelVariantExtractionStage()

        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tmp:
            pass  # Empty file

            chunk_path = Path(tmp.name)
            assert stage._validate_existing_chunk(chunk_path) is False

            # Cleanup
            chunk_path.unlink()

    def test_validate_existing_chunk_nonexistent(self):
        """Test validation rejects nonexistent files."""
        stage = ParallelVariantExtractionStage()
        chunk_path = Path("/tmp/nonexistent_file.vcf.gz")
        assert stage._validate_existing_chunk(chunk_path) is False

    def test_validate_existing_chunk_invalid_header(self):
        """Test validation rejects files with invalid VCF headers."""
        stage = ParallelVariantExtractionStage()

        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as tmp:
            # Create file with invalid header
            with gzip.open(tmp.name, "wt") as f:
                f.write("This is not a VCF file\n")
                f.write("Invalid content\n")

            chunk_path = Path(tmp.name)
            assert stage._validate_existing_chunk(chunk_path) is False

            # Cleanup
            chunk_path.unlink()

    def test_handle_checkpoint_skip(self):
        """Test checkpoint skip handler restores context correctly."""
        stage = ParallelVariantExtractionStage()

        # Mock context and workspace
        context = Mock(spec=PipelineContext)
        context.workspace = Mock()
        context.workspace.base_name = "test"
        context.workspace.get_intermediate_path = Mock()

        expected_vcf = Path("/tmp/test.variants.vcf.gz")
        context.workspace.get_intermediate_path.return_value = expected_vcf

        # Mock the file exists
        with patch.object(Path, "exists", return_value=True):
            result = stage._handle_checkpoint_skip(context)

        # Verify context was updated
        assert result == context
        assert context.extracted_vcf == expected_vcf
        assert context.data == expected_vcf
        context.mark_complete.assert_called_once_with("variant_extraction")


class TestParallelCompleteProcessingStage:
    """Test parallel complete processing stage resume functionality."""

    def test_validate_chunk_tsv_valid_file(self):
        """Test validation of valid TSV chunk."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False) as tmp:
            # Create a valid TSV file
            with gzip.open(tmp.name, "wt") as f:
                f.write("CHROM\tPOS\tREF\tALT\tSAMPLE1\n")
                f.write("chr1\t1000\tA\tT\thet\n")
                f.write("chr1\t2000\tG\tC\tref\n")

            chunk_path = Path(tmp.name)
            assert stage._validate_chunk_tsv(chunk_path) is True

            # Cleanup
            chunk_path.unlink()

    def test_validate_chunk_tsv_header_only(self):
        """Test validation rejects files with only header."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.NamedTemporaryFile(suffix=".tsv.gz", delete=False) as tmp:
            # Create file with only header
            with gzip.open(tmp.name, "wt") as f:
                f.write("CHROM\tPOS\tREF\tALT\n")

            chunk_path = Path(tmp.name)
            assert stage._validate_chunk_tsv(chunk_path) is False

            # Cleanup
            chunk_path.unlink()

    def test_validate_chunk_tsv_invalid_format(self):
        """Test validation rejects files without tab separators."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
            # Create file without tabs
            with open(tmp.name, "w") as f:
                f.write("This is not a TSV file\n")
                f.write("No tabs here either\n")

            chunk_path = Path(tmp.name)
            assert stage._validate_chunk_tsv(chunk_path) is False

            # Cleanup
            chunk_path.unlink()

    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage"
        "._validate_chunk_tsv"
    )
    def test_handle_checkpoint_skip_with_validation(self, mock_validate):
        """Test checkpoint skip handler with TSV validation."""
        stage = ParallelCompleteProcessingStage()
        mock_validate.return_value = True

        # Mock context and workspace
        context = Mock(spec=PipelineContext)
        context.workspace = Mock()
        context.workspace.base_name = "test"
        context.workspace.intermediate_dir = Path("/tmp")

        expected_tsv = Path("/tmp/test.extracted.tsv.gz")

        # Mock file exists
        with patch.object(Path, "exists", return_value=True):
            result = stage._handle_checkpoint_skip(context)

        # Verify context was updated
        assert result == context
        assert context.extracted_tsv == expected_tsv
        assert context.data == expected_tsv

        # Verify all stages marked complete
        expected_calls = [
            ("variant_extraction",),
            ("parallel_variant_extraction",),
            ("snpsift_filtering",),
            ("field_extraction",),
        ]
        actual_calls = [call.args for call in context.mark_complete.call_args_list]
        assert actual_calls == expected_calls

    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage"
        "._validate_chunk_tsv"
    )
    @patch("variantcentrifuge.stages.processing_stages.ProcessPoolExecutor")
    def test_process_chunks_parallel_with_existing_chunks(self, mock_executor_class, mock_validate):
        """Test parallel processing skips existing valid chunks."""
        stage = ParallelCompleteProcessingStage()

        # Mock context
        context = Mock(spec=PipelineContext)
        context.config = {
            "vcf_file": "/tmp/test.vcf.gz",
            "threads": 4,
            "gzip_intermediates": True,
            "fields_to_extract": "CHROM POS REF ALT",
        }
        context.workspace = Mock()
        context.workspace.base_name = "test"
        context.workspace.intermediate_dir = Path("/tmp")

        # Mock bed chunks
        bed_chunks = [Path("/tmp/chunk_0.bed"), Path("/tmp/chunk_1.bed")]

        # Mock that first chunk exists, second doesn't
        def mock_validate_side_effect(path):
            return "chunk_0" in str(path)

        mock_validate.side_effect = mock_validate_side_effect

        # Mock the executor to avoid multiprocessing issues in tests
        mock_executor = Mock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor

        # Mock future and result
        mock_future = Mock()
        mock_future.result.return_value = Path("/tmp/test.chunk_1.extracted.tsv.gz")
        mock_executor.submit.return_value = mock_future

        # Mock as_completed to return the future
        with patch(
            "variantcentrifuge.stages.processing_stages.as_completed", return_value=[mock_future]
        ):
            result = stage._process_chunks_parallel(context, bed_chunks)

        # Verify only missing chunk was submitted for processing
        assert mock_executor.submit.call_count == 1

        # Verify result contains both chunks in order
        expected_chunks = [
            Path("/tmp/test.chunk_0.extracted.tsv.gz"),
            Path("/tmp/test.chunk_1.extracted.tsv.gz"),
        ]
        assert result == expected_chunks

    def test_process_chunks_parallel_all_chunks_exist(self):
        """Test parallel processing when all chunks already exist."""
        stage = ParallelCompleteProcessingStage()

        # Mock context
        context = Mock(spec=PipelineContext)
        context.config = {
            "vcf_file": "/tmp/test.vcf.gz",
            "threads": 4,
            "gzip_intermediates": True,
            "fields_to_extract": "CHROM POS REF ALT",
        }
        context.workspace = Mock()
        context.workspace.base_name = "test"
        context.workspace.intermediate_dir = Path("/tmp")

        # Mock bed chunks
        bed_chunks = [Path("/tmp/chunk_0.bed"), Path("/tmp/chunk_1.bed")]

        # Mock that all chunks exist
        with patch.object(stage, "_validate_chunk_tsv", return_value=True):
            with patch.object(stage, "_process_single_chunk") as mock_process:
                result = stage._process_chunks_parallel(context, bed_chunks)

        # Verify no chunks were processed
        assert mock_process.call_count == 0

        # Verify result contains all chunks
        expected_chunks = [
            Path("/tmp/test.chunk_0.extracted.tsv.gz"),
            Path("/tmp/test.chunk_1.extracted.tsv.gz"),
        ]
        assert result == expected_chunks


if __name__ == "__main__":
    pytest.main([__file__])

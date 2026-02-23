"""Unit tests for ParallelCompleteProcessingStage."""

import tempfile
from concurrent.futures import Future
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest

from tests.mocks.fixtures import create_test_context
from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage


class TestParallelCompleteProcessingStage:
    """Test ParallelCompleteProcessingStage."""

    @pytest.fixture
    def context(self):
        """Create test context with config for parallel processing."""
        ctx = create_test_context()
        ctx.config["threads"] = 4
        ctx.config["vcf_file"] = "/path/to/test.vcf"
        ctx.config["extract"] = ["CHROM", "POS", "REF", "ALT"]
        ctx.config["filters"] = "QUAL >= 30"
        ctx.gene_bed_file = Path("/tmp/genes.bed")
        ctx.workspace.base_name = "test_output"
        ctx.workspace.intermediate_dir = Path("/tmp/intermediate")
        ctx.mark_complete("configuration_loading")
        ctx.mark_complete("gene_bed_creation")
        return ctx

    def test_single_thread_fallback(self, context):
        """Test that single thread falls back to sequential processing."""
        context.config["threads"] = 1
        stage = ParallelCompleteProcessingStage()

        result = stage(context)

        # Should return context unchanged for single thread
        assert result is context
        # Should not mark stages as complete
        assert not context.is_complete("variant_extraction")

    def test_dependencies(self):
        """Test stage dependencies."""
        stage = ParallelCompleteProcessingStage()
        assert stage.dependencies == {"gene_bed_creation", "configuration_loading"}
        assert stage.parallel_safe is False

    @patch("variantcentrifuge.stages.processing_stages.split_bed_file")
    @patch("variantcentrifuge.stages.processing_stages.ProcessPoolExecutor")
    @patch("variantcentrifuge.stages.processing_stages.extract_variants")
    @patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
    @patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
    def test_parallel_processing(
        self, mock_extract_fields, mock_filter, mock_extract, mock_executor, mock_split, context
    ):
        """Test parallel processing with multiple chunks."""
        # Mock BED file splitting
        chunk_paths = [Path(f"/tmp/chunks/chunk_{i}.bed") for i in range(4)]
        mock_split.return_value = [str(p) for p in chunk_paths]

        # Mock ProcessPoolExecutor
        mock_pool = MagicMock()
        mock_executor.return_value.__enter__.return_value = mock_pool

        # Create mock futures
        futures = []
        for i in range(4):
            future = Mock(spec=Future)
            future.result.return_value = Path(
                f"/tmp/intermediate/test_output.chunk_{i}.extracted.tsv.gz"
            )
            futures.append(future)

        mock_pool.submit.side_effect = futures

        # Mock as_completed to return futures in order
        with patch("variantcentrifuge.stages.processing_stages.as_completed", return_value=futures):
            stage = ParallelCompleteProcessingStage()

            # Create temp directories for test
            with tempfile.TemporaryDirectory() as tmpdir:
                tmppath = Path(tmpdir)
                chunks_dir = tmppath / "chunks"
                chunks_dir.mkdir()

                # Create mock chunk files
                for chunk_path in chunk_paths:
                    chunks_dir.joinpath(chunk_path.name).touch()

                # Create mock TSV files
                for i in range(4):
                    tsv_path = tmppath / f"test_output.chunk_{i}.extracted.tsv.gz"
                    tsv_path.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tG\n")

                # Patch cleanup to avoid deleting our test files
                with (
                    patch.object(stage, "_cleanup_chunks"),
                    patch.object(stage, "_merge_tsv_outputs", return_value=tmppath / "merged.tsv"),
                ):
                    result = stage(context)

        # Verify results
        assert result is context
        assert context.is_complete("variant_extraction")
        assert context.is_complete("snpsift_filtering")
        assert context.is_complete("field_extraction")

        # Verify split_bed_file was called
        mock_split.assert_called_once()
        assert mock_split.call_args[0][1] == 4  # n_chunks = threads

        # Verify pool was used
        assert mock_pool.submit.call_count == 4

    def test_process_single_chunk(self, context):
        """Test processing a single chunk."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            chunk_bed = tmppath / "chunk_0.bed"
            chunk_bed.write_text("chr1\t100\t200\tGENE1\n")

            config = {
                "threads_per_chunk": 2,
                "bcftools_prefilter": "QUAL > 20",
                "late_filtering": False,
                "filters": "QUAL >= 30",
                "fields_to_extract": "CHROM POS REF ALT",
                "extract_fields_separator": ",",
            }

            # Create expected output files
            chunk_vcf = tmppath / "test_output.chunk_0.variants.vcf.gz"
            chunk_vcf.touch()
            filtered_vcf = tmppath / "test_output.chunk_0.filtered.vcf.gz"
            filtered_vcf.touch()
            chunk_tsv = tmppath / "test_output.chunk_0.extracted.tsv.gz"
            chunk_tsv.touch()

            with (
                patch("variantcentrifuge.stages.processing_stages.extract_variants"),
                patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter"),
                patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools"),
            ):
                result = stage._process_single_chunk(
                    0, chunk_bed, "/path/to/test.vcf", "test_output", tmppath, config
                )

            assert result == chunk_tsv

    def test_merge_tsv_outputs_single_file(self, context):
        """Test merging when there's only one TSV file."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            single_tsv = tmppath / "chunk_0.tsv.gz"
            single_tsv.write_text("CHROM\tPOS\nchr1\t100\n")

            # Create intermediate directory
            intermediate_dir = tmppath / "intermediate"
            intermediate_dir.mkdir()

            # Mock workspace paths
            context.workspace.intermediate_dir = intermediate_dir
            context.workspace.get_intermediate_path = lambda name: intermediate_dir / name

            # Mock shutil.move
            with patch("shutil.move") as mock_move:
                stage._merge_tsv_outputs(context, [single_tsv])

            # Check that move was called with correct paths
            # The output should preserve the .gz extension from the input
            expected_output = intermediate_dir / "test_output.extracted.tsv.gz"
            mock_move.assert_called_once_with(str(single_tsv), str(expected_output))

    def test_merge_tsv_outputs_multiple_files(self, context):
        """Test merging multiple TSV files."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)

            # Create test TSV files
            tsv1 = tmppath / "chunk_0.tsv"
            tsv1.write_text("CHROM\tPOS\tREF\tALT\nchr1\t100\tA\tG\n")

            tsv2 = tmppath / "chunk_1.tsv"
            tsv2.write_text("CHROM\tPOS\tREF\tALT\nchr1\t200\tC\tT\n")

            # Create intermediate directory
            intermediate_dir = tmppath / "intermediate"
            intermediate_dir.mkdir()

            # Mock workspace paths
            context.workspace.intermediate_dir = intermediate_dir
            context.workspace.get_intermediate_path = lambda name: intermediate_dir / name
            context.config["gzip_intermediates"] = False

            result = stage._merge_tsv_outputs(context, [tsv1, tsv2])

            # Check merged content
            merged_content = result.read_text()
            lines = merged_content.strip().split("\n")
            assert len(lines) == 3  # Header + 2 data lines
            assert lines[0] == "CHROM\tPOS\tREF\tALT"
            assert "chr1\t100\tA\tG" in merged_content
            assert "chr1\t200\tC\tT" in merged_content

    def test_cleanup_chunks(self, context):
        """Test cleanup of temporary files."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            chunks_dir = tmppath / "chunks"
            chunks_dir.mkdir()

            # Create test files
            bed_chunks = [chunks_dir / f"chunk_{i}.bed" for i in range(2)]
            tsv_chunks = [tmppath / f"chunk_{i}.tsv" for i in range(2)]

            for chunk in bed_chunks + tsv_chunks:
                chunk.touch()

            # Test cleanup when keep_intermediates is False
            context.config["keep_intermediates"] = False
            stage._cleanup_chunks(bed_chunks, tsv_chunks, context)

            # All files should be deleted
            for chunk in bed_chunks + tsv_chunks:
                assert not chunk.exists()

            # Chunks directory should be removed if empty
            assert not chunks_dir.exists()

    def test_cleanup_chunks_keep_intermediates(self, context):
        """Test that cleanup preserves files when keep_intermediates is True."""
        stage = ParallelCompleteProcessingStage()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)

            # Create test files
            bed_chunk = tmppath / "chunk_0.bed"
            tsv_chunk = tmppath / "chunk_0.tsv"
            bed_chunk.touch()
            tsv_chunk.touch()

            # Test cleanup when keep_intermediates is True
            context.config["keep_intermediates"] = True
            stage._cleanup_chunks([bed_chunk], [tsv_chunk], context)

            # Files should still exist
            assert bed_chunk.exists()
            assert tsv_chunk.exists()

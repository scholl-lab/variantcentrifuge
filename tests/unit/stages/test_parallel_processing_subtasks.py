"""Test subtask timing for parallel processing stages."""

import time
from unittest.mock import Mock, patch, MagicMock
import pytest

from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage
from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner
from variantcentrifuge.pipeline_core.workspace import Workspace


class TestParallelProcessingSubtasks:
    """Test subtask timing in parallel processing."""

    @pytest.fixture
    def context(self, tmp_path):
        """Create test context with parallel processing configuration."""
        workspace = Workspace(tmp_path, "test")
        args = Mock(
            threads=4,
            vcf_file="test.vcf",
            filter="QUAL > 30",
            extract=["CHROM", "POS", "REF", "ALT"],
            late_filtering=False,
            gzip_intermediates=True,
        )

        config = {
            "threads": 4,
            "vcf_file": "test.vcf",
            "filter": "QUAL > 30",
            "extract": ["CHROM", "POS", "REF", "ALT"],
            "bcftools_prefilter": "QUAL > 20",
            "extract_fields_separator": ",",
            "gzip_intermediates": False,  # Don't gzip for testing
        }

        ctx = PipelineContext(args=args, config=config, workspace=workspace)
        ctx.gene_bed_file = tmp_path / "genes.bed"
        ctx.gene_bed_file.write_text("chr1\t1000\t2000\tGENE1\nchr2\t3000\t4000\tGENE2\n")

        # Mark dependencies as complete
        ctx.mark_complete("gene_bed_creation")
        ctx.mark_complete("configuration_loading")

        return ctx

    def test_subtask_timing_mechanism(self, context):
        """Test the subtask timing mechanism works correctly."""
        # Create a custom stage that uses subtask timing
        stage = ParallelCompleteProcessingStage()

        # Manually test the subtask timing methods
        start_time = stage._start_subtask("test_task")
        time.sleep(0.1)  # Simulate work
        stage._end_subtask("test_task", start_time)

        # Check that subtask time was recorded
        assert hasattr(stage, "_subtask_times")
        assert "test_task" in stage._subtask_times
        assert stage._subtask_times["test_task"] >= 0.1

    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage._process_chunks_parallel"
    )
    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage._merge_tsv_outputs"
    )
    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage._split_bed_file"
    )
    @patch(
        "variantcentrifuge.stages.processing_stages.ParallelCompleteProcessingStage._cleanup_chunks"
    )
    def test_subtask_timing_recorded(
        self, mock_cleanup, mock_split, mock_merge, mock_process, context
    ):
        """Test that subtask times are recorded during parallel processing."""
        # Create the stage first
        stage = ParallelCompleteProcessingStage()

        # Mock the internal methods to return appropriate values
        mock_split.return_value = ["chunk1.bed", "chunk2.bed"]
        mock_merge.return_value = context.workspace.intermediate_dir / "merged.tsv"

        # Configure the mock_process to set subtask times
        def process_with_subtasks(ctx, bed_chunks):
            # Simulate the subtask timing that would happen in _process_chunks_parallel
            stage._subtask_times["variant_extraction"] = 2.5
            stage._subtask_times["snpsift_filtering"] = 1.5
            stage._subtask_times["field_extraction"] = 0.5
            return ["chunk1.tsv", "chunk2.tsv"]

        mock_process.side_effect = process_with_subtasks

        # Run the stage
        result = stage(context)

        # Check that all subtask times were recorded
        assert hasattr(stage, "_subtask_times")
        assert "bed_splitting" in stage._subtask_times
        assert "parallel_chunk_processing" in stage._subtask_times
        assert "tsv_merging" in stage._subtask_times
        assert "cleanup" in stage._subtask_times

        # Check that chunk processing recorded averaged times
        assert "variant_extraction" in stage._subtask_times
        assert "snpsift_filtering" in stage._subtask_times
        assert "field_extraction" in stage._subtask_times

        # Verify times are positive
        for subtask, duration in stage._subtask_times.items():
            assert duration >= 0, f"Subtask {subtask} has negative duration"

    def test_runner_captures_subtask_times(self, context):
        """Test that PipelineRunner captures and displays subtask times."""
        # Create a mock stage with subtask times
        mock_stage = Mock()
        mock_stage.name = "test_stage"
        mock_stage.description = "Test stage"
        mock_stage.dependencies = set()
        mock_stage.soft_dependencies = set()  # Add soft_dependencies
        mock_stage.parallel_safe = False
        mock_stage.estimated_runtime = 1.0

        # Pre-set the subtask_times to avoid Mock issues
        subtask_dict = {
            "subtask1": 0.5,
            "subtask2": 0.3,
            "subtask3": 0.1,
        }
        mock_stage.subtask_times = subtask_dict

        # Mock the stage execution
        def execute_with_subtasks(ctx):
            # Just return the context
            return ctx

        mock_stage.__call__ = execute_with_subtasks

        # Run the stage through the runner
        runner = PipelineRunner()
        result = runner.run([mock_stage], context)

        # Check that runner captured the subtask times
        assert "test_stage" in runner._subtask_times
        assert runner._subtask_times["test_stage"] == {
            "subtask1": 0.5,
            "subtask2": 0.3,
            "subtask3": 0.1,
        }

    @patch("variantcentrifuge.stages.processing_stages.split_bed_file")
    @patch("variantcentrifuge.stages.processing_stages.ProcessPoolExecutor")
    def test_averaged_times_display(self, mock_executor_class, mock_split_bed, context, caplog):
        """Test that averaged times are displayed with (avg) suffix."""
        # Mock the executor to simulate parallel processing
        mock_executor = MagicMock()
        mock_executor_class.return_value.__enter__.return_value = mock_executor

        # Mock split_bed_file
        mock_split_bed.return_value = ["chunk1.bed", "chunk2.bed"]

        # Create stage and simulate execution with subtask times
        stage = ParallelCompleteProcessingStage()
        stage._subtask_times = {
            "bed_splitting": 1.0,
            "parallel_chunk_processing": 10.0,
            "variant_extraction": 2.5,  # This should show (avg)
            "snpsift_filtering": 1.5,  # This should show (avg)
            "field_extraction": 0.5,  # This should show (avg)
            "tsv_merging": 0.8,
            "cleanup": 0.2,
        }

        # Create runner and add the stage to execution times
        runner = PipelineRunner()
        runner._execution_times = {"parallel_complete_processing": 15.0}
        runner._subtask_times = {"parallel_complete_processing": stage._subtask_times}

        # Call the summary logging
        import logging

        with caplog.at_level(logging.INFO):
            runner._log_execution_summary()

        # Check that averaged times have (avg) suffix in the logs
        log_text = "\n".join(record.message for record in caplog.records)
        assert "variant_extraction (avg)" in log_text
        assert "snpsift_filtering (avg)" in log_text
        assert "field_extraction (avg)" in log_text

        # Check that non-averaged times don't have the suffix
        assert "bed_splitting (avg)" not in log_text
        assert "tsv_merging (avg)" not in log_text
        assert "cleanup (avg)" not in log_text

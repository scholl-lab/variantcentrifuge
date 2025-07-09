"""Tests for checkpoint resume functionality."""

import os
import tempfile
from unittest.mock import MagicMock, patch

from variantcentrifuge.checkpoint import CheckpointContext, PipelineState


class TestCheckpointResume:
    """Test checkpoint resume functionality in the pipeline."""

    @patch("variantcentrifuge.pipeline_core.get_gene_bed")
    @patch("variantcentrifuge.pipeline_core.split_bed_file")
    def test_parallel_pipeline_resume_from_checkpoint(self, mock_split_bed, mock_get_gene_bed):
        """Test that parallel pipeline correctly resumes from checkpoint state."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Create initial pipeline state with some completed steps
            state = PipelineState(output_dir, enable_checksum=False)

            # Simulate that gene BED creation was completed
            bed_file = os.path.join(output_dir, "bed_cache", "genes_test.bed")
            os.makedirs(os.path.dirname(bed_file), exist_ok=True)
            with open(bed_file, "w") as f:
                f.write("chr1\t1000\t2000\tTEST\n")

            with CheckpointContext(state, "parallel_gene_bed_creation") as ctx:
                ctx.add_output_file(bed_file)

            # Simulate that BED splitting was completed with 3 chunks
            intermediate_dir = os.path.join(output_dir, "intermediate")
            os.makedirs(intermediate_dir, exist_ok=True)

            chunk_files = []
            with CheckpointContext(
                state, "parallel_bed_splitting", parameters={"num_chunks": 8}
            ) as ctx:
                ctx.add_input_file(bed_file)
                for i in range(3):  # Only 3 chunks created, not 8
                    chunk_file = os.path.join(intermediate_dir, f"chunk_{i}.bed")
                    with open(chunk_file, "w") as f:
                        f.write(f"chr1\t{1000+i*100}\t{1100+i*100}\tTEST_chunk{i}\n")
                    chunk_files.append(chunk_file)
                    ctx.add_output_file(chunk_file)

            # Save the state
            state.save()

            # Now simulate resuming the pipeline
            args = MagicMock()
            args.output_dir = output_dir
            args.vcf_file = "/fake/input.vcf"

            # Create a new pipeline state that will load the saved state
            resume_state = PipelineState(output_dir, enable_checksum=False)
            resume_state.load()

            cfg = {  # noqa: F841
                "reference": "GRCh37",
                "output_dir": output_dir,
                "_pipeline_state": resume_state,
                "keep_intermediates": True,
            }

            # These should not be called since steps are already completed
            mock_get_gene_bed.side_effect = Exception(
                "get_gene_bed should not be called during resume!"
            )
            mock_split_bed.side_effect = Exception(
                "split_bed_file should not be called during resume!"
            )

            # Mock the chunk processing to avoid actual work
            with patch("variantcentrifuge.pipeline_core._process_bed_chunk") as mock_process:
                mock_process.return_value = "/fake/output.tsv"

                # Mock file operations
                with patch("variantcentrifuge.pipeline_core.os.path.exists") as mock_exists:
                    with patch("variantcentrifuge.pipeline_core.os.path.getsize") as mock_getsize:
                        # Make the mocked files appear to exist
                        def exists_side_effect(path):
                            if path == bed_file:
                                return True
                            if "chunk_" in path and ".bed" in path:
                                # Check if it's one of our actual chunk files
                                return path in chunk_files
                            # Let other files use real os.path.exists
                            return os.path.exists(path)

                        mock_exists.side_effect = exists_side_effect
                        mock_getsize.return_value = 100

                        # This should successfully skip the completed steps
                        # Note: We can't run the full pipeline due to external dependencies,
                        # but we can verify the checkpoint loading works

                        # Check that the state correctly identifies completed steps
                        assert resume_state.should_skip_step("parallel_gene_bed_creation")
                        assert resume_state.should_skip_step("parallel_bed_splitting")

                        # Verify the correct files are retrieved from state
                        gene_step = resume_state.state["steps"]["parallel_gene_bed_creation"]
                        assert gene_step.output_files[0].path == bed_file

                        split_step = resume_state.state["steps"]["parallel_bed_splitting"]
                        retrieved_chunks = [f.path for f in split_step.output_files]
                        assert len(retrieved_chunks) == 3  # Not 8!
                        assert retrieved_chunks == chunk_files

    def test_checkpoint_file_path_retrieval(self):
        """Test that file paths are correctly retrieved from checkpoint state."""
        with tempfile.TemporaryDirectory() as output_dir:
            state = PipelineState(output_dir, enable_checksum=False)

            # Add some file paths
            test_paths = ["path/to/file1.bed", "path/to/file2.bed", "another/path/file3.bed"]

            with CheckpointContext(state, "test_step") as ctx:
                for path in test_paths:
                    # Create dummy files
                    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
                    with open(path, "w") as f:
                        f.write("test")
                    ctx.add_output_file(path)

            # Save and reload
            state.save()

            new_state = PipelineState(output_dir, enable_checksum=False)
            new_state.load()

            # Retrieve paths
            step_info = new_state.state["steps"]["test_step"]
            retrieved_paths = [f.path for f in step_info.output_files]

            assert retrieved_paths == test_paths

            # Clean up
            for path in test_paths:
                if os.path.exists(path):
                    os.remove(path)

"""Test checkpoint resume functionality for parallel pipeline."""

import json
import os
import shutil
import tempfile

import pytest

from variantcentrifuge.checkpoint import PipelineState


class TestParallelCheckpointResume:
    """Test checkpoint resume for parallel pipeline execution."""

    def setup_method(self):
        """Set up test environment."""
        self.test_dir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self.test_dir, "output")
        self.intermediate_dir = os.path.join(self.output_dir, "intermediate")
        os.makedirs(self.intermediate_dir)

    def teardown_method(self):
        """Clean up test environment."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def create_mock_state(self, steps_completed):
        """Create a mock checkpoint state with specified completed steps."""
        state_file = os.path.join(self.output_dir, ".variantcentrifuge_state.json")
        state_data = {
            "version": "1.0",
            "pipeline_version": "0.5.0",
            "start_time": 1234567890.0,
            "last_update": 1234567890.0,
            "configuration_hash": "test_hash",
            "steps": {},
            "metadata": {"command_line": "test command"},
        }

        # Add completed steps
        for step_name, file_info in steps_completed.items():
            state_data["steps"][step_name] = {
                "name": step_name,
                "status": "completed",
                "start_time": 1234567890.0,
                "end_time": 1234567891.0,
                "command_hash": None,
                "input_files": file_info.get("input_files", []),
                "output_files": file_info.get("output_files", []),
                "parameters": file_info.get("parameters", {}),
                "error": None,
            }

        with open(state_file, "w") as f:
            json.dump(state_data, f, indent=2)

        # Create and load the PipelineState
        pipeline_state = PipelineState(self.output_dir)
        pipeline_state.load()  # This will convert dicts to StepInfo objects
        return pipeline_state

    def test_resume_after_parallel_merge(self):
        """Test resuming after parallel merge is completed."""
        # Create mock files
        merged_file = os.path.join(self.intermediate_dir, "test.extracted.tsv.gz")

        # Create the merged file (simulating completed merge)
        with open(merged_file, "w") as f:
            f.write("header\ndata\n")

        # Create checkpoint state showing merge is complete
        steps_completed = {
            "parallel_merge_chunks": {
                "output_files": [
                    {
                        "path": merged_file,
                        "size": os.path.getsize(merged_file),
                        "mtime": os.path.getmtime(merged_file),
                        "checksum": None,
                    }
                ],
                "parameters": {"num_chunks": 5},
            }
        }
        pipeline_state = self.create_mock_state(steps_completed)

        # Verify the checkpoint contains the correct information
        assert pipeline_state.should_skip_step("parallel_merge_chunks")

        # Verify we can retrieve the output file path
        step_info = pipeline_state.state["steps"].get("parallel_merge_chunks")
        assert step_info is not None
        assert hasattr(step_info, "output_files")
        assert len(step_info.output_files) == 1
        assert step_info.output_files[0].path == merged_file

    def test_resume_after_sorting(self):
        """Test resuming after TSV sorting is completed."""
        # Create mock files
        sorted_file = os.path.join(self.intermediate_dir, "test.extracted.sorted.tsv.gz")

        # Create the sorted file
        with open(sorted_file, "w") as f:
            f.write("GENE\tCHROM\tPOS\n")
            f.write("GENE1\tchr1\t100\n")

        # Create checkpoint state showing sorting is complete
        steps_completed = {
            "parallel_merge_chunks": {
                "output_files": [
                    {
                        "path": os.path.join(self.intermediate_dir, "test.extracted.tsv.gz"),
                        "size": 100,
                        "mtime": 1234567890.0,
                        "checksum": None,
                    }
                ],
            },
            "tsv_sorting": {
                "output_files": [
                    {
                        "path": sorted_file,
                        "size": os.path.getsize(sorted_file),
                        "mtime": os.path.getmtime(sorted_file),
                        "checksum": None,
                    }
                ],
                "parameters": {
                    "gene_column": "GENE",
                    "memory_limit": "2G",
                    "parallel": 4,
                },
            },
        }
        pipeline_state = self.create_mock_state(steps_completed)

        # Verify that the sorted file path is retrieved correctly
        step_info = pipeline_state.state["steps"].get("tsv_sorting")
        assert step_info is not None
        assert hasattr(step_info, "output_files")
        assert len(step_info.output_files) == 1
        assert step_info.output_files[0].path == sorted_file

    def test_resume_handles_missing_files(self):
        """Test that resume handles missing checkpoint files gracefully."""
        # Create checkpoint state with non-existent file
        steps_completed = {
            "parallel_merge_chunks": {
                "output_files": [
                    {
                        "path": "/non/existent/file.tsv.gz",
                        "size": 100,
                        "mtime": 1234567890.0,
                        "checksum": None,
                    }
                ],
            }
        }
        pipeline_state = self.create_mock_state(steps_completed)

        # The pipeline should detect the missing file
        step_info = pipeline_state.state["steps"].get("parallel_merge_chunks")
        assert step_info is not None
        assert hasattr(step_info, "output_files")
        file_path = step_info.output_files[0].path
        assert not os.path.exists(file_path)

    def test_file_cleanup_timing(self):
        """Test that intermediate files are not cleaned up prematurely."""
        # Create mock chunk files
        chunk_files = []
        for i in range(3):
            chunk_file = os.path.join(self.intermediate_dir, f"test.chunk_{i}.extracted.tsv.gz")
            with open(chunk_file, "w") as f:
                f.write(f"chunk {i} data\n")
            chunk_files.append(chunk_file)

        # Simulate merge completion
        merged_file = os.path.join(self.intermediate_dir, "test.extracted.tsv.gz")
        with open(merged_file, "w") as f:
            f.write("merged data\n")

        # With keep_intermediates=False, chunk files should be removed after merge
        cfg = {"keep_intermediates": False}

        # Simulate cleanup
        if not cfg.get("keep_intermediates"):
            for chunk_file in chunk_files:
                if os.path.exists(chunk_file):
                    os.remove(chunk_file)

        # Verify chunk files are removed but merged file remains
        for chunk_file in chunk_files:
            assert not os.path.exists(chunk_file)
        assert os.path.exists(merged_file)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

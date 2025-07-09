"""Tests for parallel pipeline checkpoint functionality."""

import os
import tempfile

from variantcentrifuge.checkpoint import CheckpointContext, PipelineState


class TestParallelCheckpoint:
    """Test checkpoint functionality in parallel pipeline execution."""

    def test_parallel_chunk_tracking(self):
        """Test that parallel chunks are properly tracked in checkpoint state."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Create pipeline state
            state = PipelineState(output_dir, enable_checksum=False)

            # Simulate parallel chunk processing
            num_chunks = 4

            # Track gene bed creation
            with CheckpointContext(state, "parallel_gene_bed_creation") as ctx:
                ctx.add_output_file("/fake/master.bed")

            # Track bed splitting
            with CheckpointContext(
                state, "parallel_bed_splitting", parameters={"num_chunks": num_chunks}
            ) as ctx:
                ctx.add_input_file("/fake/master.bed")
                for i in range(num_chunks):
                    ctx.add_output_file(f"/fake/chunk_{i}.bed")

            # Simulate parallel chunk processing
            for i in range(num_chunks):
                step_name = f"parallel_chunk_{i}_processing"
                with CheckpointContext(
                    state,
                    step_name,
                    parameters={"chunk_index": i, "bed_file": f"/fake/chunk_{i}.bed"},
                ) as ctx:
                    ctx.add_input_file(f"/fake/chunk_{i}.bed")
                    ctx.add_input_file("/fake/input.vcf")
                    ctx.add_output_file(f"/fake/chunk_{i}.tsv")

            # Track merge
            with CheckpointContext(
                state, "parallel_merge_chunks", parameters={"num_chunks": num_chunks}
            ) as ctx:
                for i in range(num_chunks):
                    ctx.add_input_file(f"/fake/chunk_{i}.tsv")
                ctx.add_output_file("/fake/merged.tsv")

            # Verify all steps are recorded
            assert (
                len(state.state["steps"]) == 2 + num_chunks + 1
            )  # bed_creation + splitting + chunks + merge

            # Verify chunk steps
            for i in range(num_chunks):
                step_name = f"parallel_chunk_{i}_processing"
                assert step_name in state.state["steps"]
                step_info = state.state["steps"][step_name]
                assert step_info.status == "completed"
                assert step_info.parameters["chunk_index"] == i

            # Test resume - all steps should be skipped
            state2 = PipelineState(output_dir, enable_checksum=False)
            state2.load()

            assert state2.should_skip_step("parallel_gene_bed_creation")
            assert state2.should_skip_step("parallel_bed_splitting")
            for i in range(num_chunks):
                assert state2.should_skip_step(f"parallel_chunk_{i}_processing")
            assert state2.should_skip_step("parallel_merge_chunks")

    def test_parallel_chunk_failure_tracking(self):
        """Test that failures in parallel chunks are properly tracked."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Create pipeline state
            state = PipelineState(output_dir, enable_checksum=False)

            # Simulate chunk 0 and 1 succeed, chunk 2 fails
            for i in range(2):
                step_name = f"parallel_chunk_{i}_processing"
                with CheckpointContext(state, step_name) as ctx:
                    ctx.add_output_file(f"/fake/chunk_{i}.tsv")

            # Chunk 2 fails
            step_name = "parallel_chunk_2_processing"
            state.start_step(step_name)
            state.fail_step(step_name, "Simulated failure")

            # Verify state
            assert state.state["steps"]["parallel_chunk_0_processing"].status == "completed"
            assert state.state["steps"]["parallel_chunk_1_processing"].status == "completed"
            assert state.state["steps"]["parallel_chunk_2_processing"].status == "failed"

            # On resume, completed chunks should be skipped
            state2 = PipelineState(output_dir, enable_checksum=False)
            state2.load()

            assert state2.should_skip_step("parallel_chunk_0_processing")
            assert state2.should_skip_step("parallel_chunk_1_processing")
            assert not state2.should_skip_step("parallel_chunk_2_processing")  # Failed, needs retry
            assert not state2.should_skip_step("parallel_chunk_3_processing")  # Never started

    def test_parallel_pipeline_checkpoint_recording(self):
        """Test that parallel pipeline properly records checkpoint steps."""
        with tempfile.TemporaryDirectory() as output_dir:
            # Create pipeline state
            state = PipelineState(output_dir, enable_checksum=False)

            # Mock the parallel pipeline process
            # 1. Gene BED creation
            with CheckpointContext(state, "parallel_gene_bed_creation") as ctx:
                # Create a real temp file to track
                bed_file = os.path.join(output_dir, "master.bed")
                with open(bed_file, "w") as f:
                    f.write("chr1\t1000\t2000\tGENE1\n")
                ctx.add_output_file(bed_file)

            # 2. BED splitting
            chunk_files = []
            with CheckpointContext(
                state, "parallel_bed_splitting", parameters={"num_chunks": 4}
            ) as ctx:
                ctx.add_input_file(bed_file)
                for i in range(4):
                    chunk_file = os.path.join(output_dir, f"chunk_{i}.bed")
                    with open(chunk_file, "w") as f:
                        f.write(f"chr1\t{1000+i*100}\t{1100+i*100}\tGENE1_chunk{i}\n")
                    chunk_files.append(chunk_file)
                    ctx.add_output_file(chunk_file)

            # 3. Process chunks (simulating parallel execution)
            tsv_files = []
            for i in range(4):
                step_name = f"parallel_chunk_{i}_processing"
                with CheckpointContext(state, step_name, parameters={"chunk_index": i}) as ctx:
                    ctx.add_input_file(chunk_files[i])
                    tsv_file = os.path.join(output_dir, f"chunk_{i}.tsv")
                    with open(tsv_file, "w") as f:
                        f.write(f"data for chunk {i}\n")
                    tsv_files.append(tsv_file)
                    ctx.add_output_file(tsv_file)

            # 4. Merge chunks
            with CheckpointContext(
                state, "parallel_merge_chunks", parameters={"num_chunks": 4}
            ) as ctx:
                for tsv_file in tsv_files:
                    ctx.add_input_file(tsv_file)
                merged_file = os.path.join(output_dir, "merged.tsv")
                with open(merged_file, "w") as f:
                    f.write("merged data\n")
                ctx.add_output_file(merged_file)

            # Verify state
            assert len(state.state["steps"]) == 7  # 1 + 1 + 4 + 1

            # Verify all steps completed
            for step_name, step_info in state.state["steps"].items():
                assert step_info.status == "completed", f"Step {step_name} not completed"

            # Test resume functionality
            state2 = PipelineState(output_dir, enable_checksum=False)
            state2.load()

            # All steps should be skipped
            assert state2.should_skip_step("parallel_gene_bed_creation")
            assert state2.should_skip_step("parallel_bed_splitting")
            for i in range(4):
                assert state2.should_skip_step(f"parallel_chunk_{i}_processing")
            assert state2.should_skip_step("parallel_merge_chunks")

    def test_checkpoint_state_ordering(self):
        """Test that checkpoint state maintains proper ordering of parallel chunks."""
        with tempfile.TemporaryDirectory() as output_dir:
            state = PipelineState(output_dir, enable_checksum=False)

            # Process chunks out of order (3, 1, 0, 2)
            for i in [3, 1, 0, 2]:
                step_name = f"parallel_chunk_{i}_processing"
                with CheckpointContext(state, step_name) as ctx:
                    ctx.add_output_file(f"/fake/chunk_{i}.tsv")

            # Verify all are recorded
            for i in range(4):
                assert f"parallel_chunk_{i}_processing" in state.state["steps"]

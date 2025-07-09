"""Tests for the checkpoint and resume system."""

import os
import time

from variantcentrifuge.checkpoint import (
    CheckpointContext,
    FileInfo,
    PipelineState,
    StepInfo,
    checkpoint,
)


class TestFileInfo:
    """Test FileInfo class."""

    def test_file_info_creation(self, tmp_path):
        """Test creating FileInfo from a file."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")

        # Without checksum
        info = FileInfo.from_file(str(test_file), calculate_checksum=False)
        assert info.path == str(test_file)
        assert info.size == 12  # "test content" is 12 bytes
        assert info.mtime > 0
        assert info.checksum is None

        # With checksum
        info_with_checksum = FileInfo.from_file(str(test_file), calculate_checksum=True)
        assert info_with_checksum.checksum is not None
        assert len(info_with_checksum.checksum) == 64  # SHA256 hex digest

    def test_file_info_validation(self, tmp_path):
        """Test FileInfo validation."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")

        info = FileInfo.from_file(str(test_file))

        # Should validate successfully
        assert info.validate() is True

        # Modify file - validation should fail
        test_file.write_text("modified content")
        assert info.validate() is False

        # Delete file - validation should fail
        test_file.unlink()
        assert info.validate() is False

    def test_file_info_checksum_validation(self, tmp_path):
        """Test FileInfo checksum validation."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")

        info = FileInfo.from_file(str(test_file), calculate_checksum=True)
        assert info.validate(calculate_checksum=True) is True

        # Modify file content
        test_file.write_text("different content")
        assert info.validate(calculate_checksum=True) is False


class TestStepInfo:
    """Test StepInfo class."""

    def test_step_info_creation(self):
        """Test creating StepInfo."""
        step = StepInfo(name="test_step", status="pending")
        assert step.name == "test_step"
        assert step.status == "pending"
        assert step.duration is None

    def test_step_info_duration(self):
        """Test step duration calculation."""
        step = StepInfo(name="test_step", status="completed", start_time=100.0, end_time=150.5)
        assert step.duration == 50.5

    def test_step_info_serialization(self, tmp_path):
        """Test StepInfo serialization to/from dict."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test")

        file_info = FileInfo.from_file(str(test_file))

        step = StepInfo(
            name="test_step",
            status="completed",
            start_time=100.0,
            end_time=150.0,
            command_hash="abc123",
            input_files=[file_info],
            output_files=[file_info],
            parameters={"param1": "value1"},
            error=None,
        )

        # Convert to dict
        step_dict = step.to_dict()
        assert step_dict["name"] == "test_step"
        assert step_dict["status"] == "completed"
        assert len(step_dict["input_files"]) == 1

        # Convert back from dict
        restored_step = StepInfo.from_dict(step_dict)
        assert restored_step.name == step.name
        assert restored_step.status == step.status
        assert restored_step.duration == step.duration
        assert len(restored_step.input_files) == 1


class TestPipelineState:
    """Test PipelineState class."""

    def test_pipeline_state_initialization(self, tmp_path):
        """Test PipelineState initialization."""
        state = PipelineState(str(tmp_path))
        assert state.output_dir == str(tmp_path)
        assert state.state_file == str(tmp_path / ".variantcentrifuge_state.json")
        assert not state._loaded_from_file

    def test_pipeline_state_save_load(self, tmp_path):
        """Test saving and loading pipeline state."""
        state = PipelineState(str(tmp_path))
        state.initialize({"test": "config"}, "1.0.0")

        # Add a step
        state.start_step("test_step", command_hash="abc123")
        state.complete_step("test_step")

        # Save state
        state.save()
        assert os.path.exists(state.state_file)

        # Load into new instance
        new_state = PipelineState(str(tmp_path))
        assert new_state.load() is True
        assert new_state._loaded_from_file is True
        assert "test_step" in new_state.state["steps"]
        assert new_state.state["steps"]["test_step"].status == "completed"

    def test_pipeline_state_resume_compatibility(self, tmp_path):
        """Test pipeline state resume compatibility checking."""
        state = PipelineState(str(tmp_path))
        # Use actual pipeline configuration keys
        config = {"gene_name": "BRCA1", "vcf_file": "input.vcf", "filters": "QUAL > 30"}
        state.initialize(config, "1.0.0")
        state.save()

        # Load and check compatibility
        new_state = PipelineState(str(tmp_path))
        new_state.load()

        # Same config and version - should be compatible
        assert new_state.can_resume(config, "1.0.0") is True

        # Different config - should not be compatible
        different_config = {"gene_name": "BRCA2", "vcf_file": "input.vcf", "filters": "QUAL > 30"}
        assert new_state.can_resume(different_config, "1.0.0") is False

        # Different version - should not be compatible
        assert new_state.can_resume(config, "2.0.0") is False

    def test_pipeline_state_step_tracking(self, tmp_path):
        """Test step tracking functionality."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        # Start a step
        state.start_step("step1", parameters={"param": "value"})
        assert state.state["steps"]["step1"].status == "running"

        # Complete the step
        test_file = tmp_path / "output.txt"
        test_file.write_text("output")
        state.complete_step("step1", output_files=[str(test_file)])

        assert state.state["steps"]["step1"].status == "completed"
        assert len(state.state["steps"]["step1"].output_files) == 1
        assert state.state["steps"]["step1"].duration > 0

        # Fail a step
        state.start_step("step2")
        state.fail_step("step2", "Test error")
        assert state.state["steps"]["step2"].status == "failed"
        assert state.state["steps"]["step2"].error == "Test error"

    def test_pipeline_state_resume_point(self, tmp_path):
        """Test finding resume point."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        # No completed steps
        assert state.get_resume_point() is None

        # Add completed steps with delays to ensure different timestamps
        state.start_step("step1")
        time.sleep(0.01)
        state.complete_step("step1")

        time.sleep(0.01)
        state.start_step("step2")
        time.sleep(0.01)
        state.complete_step("step2")

        time.sleep(0.01)
        state.start_step("step3")
        state.fail_step("step3", "Failed")

        # Should return the last completed step
        assert state.get_resume_point() == "step2"

    def test_pipeline_state_should_skip(self, tmp_path):
        """Test step skipping logic."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        # Not loaded from file - should not skip
        assert state.should_skip_step("step1") is False

        # Save and reload
        state.start_step("step1")
        state.complete_step("step1")
        state.save()

        new_state = PipelineState(str(tmp_path))
        new_state.load()

        # Completed step should be skipped
        assert new_state.should_skip_step("step1") is True

        # Non-existent step should not be skipped
        assert new_state.should_skip_step("step2") is False

    def test_pipeline_state_summary(self, tmp_path):
        """Test summary generation."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        state.start_step("step1")
        state.complete_step("step1")

        state.start_step("step2")
        state.fail_step("step2", "Test error")

        summary = state.get_summary()
        assert "Pipeline State Summary:" in summary
        assert "✓ step1" in summary
        assert "✗ step2" in summary
        assert "Error: Test error" in summary


class TestCheckpointDecorator:
    """Test checkpoint decorator."""

    def test_checkpoint_decorator_no_state(self):
        """Test checkpoint decorator without pipeline state."""

        @checkpoint("test_step")
        def test_function(value):
            return value * 2

        # Should work normally without pipeline state
        result = test_function(5)
        assert result == 10

    def test_checkpoint_decorator_with_state(self, tmp_path):
        """Test checkpoint decorator with pipeline state."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        call_count = 0

        @checkpoint("test_step", parameters={"key": "value"})
        def test_function(value):
            nonlocal call_count
            call_count += 1
            return value * 2

        # First call - should execute
        result = test_function(5, _pipeline_state=state)
        assert result == 10
        assert call_count == 1
        assert state.state["steps"]["test_step"].status == "completed"

        # Save and reload
        state.save()
        new_state = PipelineState(str(tmp_path))
        new_state.load()

        # Second call with reloaded state - should skip
        result = test_function(5, _pipeline_state=new_state)
        assert result is None  # Default return when skipped
        assert call_count == 1  # Function not called again

    def test_checkpoint_decorator_with_files(self, tmp_path):
        """Test checkpoint decorator with file tracking."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        input_file = tmp_path / "input.txt"
        output_file = tmp_path / "output.txt"
        input_file.write_text("input")

        @checkpoint("file_step", input_files=str(input_file), output_files=str(output_file))
        def process_file():
            output_file.write_text("output")
            return str(output_file)

        result = process_file(_pipeline_state=state)
        assert result == str(output_file)

        step = state.state["steps"]["file_step"]
        assert len(step.input_files) == 1
        assert len(step.output_files) == 1
        assert step.input_files[0].path == str(input_file)
        assert step.output_files[0].path == str(output_file)


class TestCheckpointContext:
    """Test CheckpointContext context manager."""

    def test_checkpoint_context_no_state(self):
        """Test checkpoint context without pipeline state."""
        with CheckpointContext(None, "test_step") as ctx:
            assert ctx.skip is False
            # Should work without errors

    def test_checkpoint_context_with_state(self, tmp_path):
        """Test checkpoint context with pipeline state."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        # First execution
        with CheckpointContext(state, "test_step", parameters={"p": 1}) as ctx:
            assert ctx.skip is False
            ctx.add_input_file(str(tmp_path / "input.txt"))
            ctx.add_output_file(str(tmp_path / "output.txt"))

        assert state.state["steps"]["test_step"].status == "completed"

        # Save and reload
        state.save()
        new_state = PipelineState(str(tmp_path))
        new_state.load()

        # Second execution - should skip
        executed = False
        with CheckpointContext(new_state, "test_step") as ctx:
            if not ctx.skip:
                executed = True

        assert executed is False

    def test_checkpoint_context_error_handling(self, tmp_path):
        """Test checkpoint context error handling."""
        state = PipelineState(str(tmp_path))
        state.initialize({}, "1.0.0")

        try:
            with CheckpointContext(state, "error_step") as ctx:
                assert ctx.skip is False
                raise ValueError("Test error")
        except ValueError:
            pass

        assert state.state["steps"]["error_step"].status == "failed"
        assert "Test error" in state.state["steps"]["error_step"].error


class TestIntegration:
    """Integration tests for checkpoint system."""

    def test_full_pipeline_simulation(self, tmp_path):
        """Test simulating a full pipeline run with checkpoints."""
        state = PipelineState(str(tmp_path), enable_checksum=True)
        state.initialize({"param": "value"}, "1.0.0")

        # Simulate pipeline steps
        steps_executed = []

        def step1():
            steps_executed.append("step1")
            output = tmp_path / "step1_output.txt"
            output.write_text("step1 result")
            return str(output)

        def step2(input_file):
            steps_executed.append("step2")
            output = tmp_path / "step2_output.txt"
            output.write_text(f"step2 processed {input_file}")
            return str(output)

        def step3(input_file):
            steps_executed.append("step3")
            # This step will fail
            raise RuntimeError("Step 3 failed")

        # Run pipeline
        with CheckpointContext(state, "step1") as ctx:
            if not ctx.skip:
                result1 = step1()
                ctx.add_output_file(result1)

        with CheckpointContext(state, "step2") as ctx:
            if not ctx.skip:
                result2 = step2(result1)
                ctx.add_input_file(result1)
                ctx.add_output_file(result2)

        try:
            with CheckpointContext(state, "step3") as ctx:
                if not ctx.skip:
                    step3(result2)
        except RuntimeError:
            pass

        assert steps_executed == ["step1", "step2", "step3"]

        # Save state
        state.save()

        # Simulate resume after fixing step3
        new_state = PipelineState(str(tmp_path), enable_checksum=True)
        new_state.load()

        steps_executed.clear()

        # Resume - should skip completed steps
        with CheckpointContext(new_state, "step1") as ctx:
            if not ctx.skip:
                result1 = step1()
            else:
                result1 = str(tmp_path / "step1_output.txt")

        with CheckpointContext(new_state, "step2") as ctx:
            if not ctx.skip:
                result2 = step2(result1)
            else:
                result2 = str(tmp_path / "step2_output.txt")

        # Fixed step3
        def step3_fixed(input_file):
            steps_executed.append("step3_fixed")
            output = tmp_path / "step3_output.txt"
            output.write_text(f"step3 processed {input_file}")
            return str(output)

        with CheckpointContext(new_state, "step3") as ctx:
            if not ctx.skip:
                result3 = step3_fixed(result2)
                ctx.add_input_file(result2)
                ctx.add_output_file(result3)

        # Only step3 should have been executed
        assert steps_executed == ["step3_fixed"]
        assert new_state.state["steps"]["step3"].status == "completed"

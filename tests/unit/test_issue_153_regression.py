"""Regression tests for Issue #153:
- PipelineRunner._should_skip_stage and Stage.enable_checkpoint
- PipelineState.should_skip_step validates recorded outputs exist and match size/mtime/checksum
- Invalidate downstream stages when upstream re-runs
- Atomic write of .variantcentrifuge_state.json
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest

from variantcentrifuge.checkpoint import PipelineState
from variantcentrifuge.pipeline_core import PipelineContext, PipelineRunner, Stage


class DummyStage(Stage):
    """Simple stage for testing checkpoint behavior."""

    def __init__(self, name: str, dependencies: list[str] | None = None, enable_chk: bool = True):
        super().__init__()
        self._name = name
        self._deps = set(dependencies or [])
        self._enable_chk = enable_chk
        self.run_count = 0

    @property
    def name(self) -> str:
        return self._name

    @property
    def dependencies(self) -> set[str]:
        return self._deps

    @property
    def enable_checkpoint(self) -> bool:
        return self._enable_chk

    @property
    def checkpoint_resume_policy(self) -> str:
        return "skip"

    def _process(self, context: PipelineContext) -> PipelineContext:
        self.run_count += 1
        out_file = Path(context.config["output_dir"]) / f"{self.name}.out"
        out_file.write_text(f"content from {self.name} run {self.run_count}\n")
        context.mark_complete(self.name)
        return context

    def get_output_files(self, context: PipelineContext) -> list[str]:
        out_file = Path(context.config["output_dir"]) / f"{self.name}.out"
        return [str(out_file)] if out_file.exists() else []


@pytest.mark.unit
def test_should_skip_step_validates_output_files_exist(tmp_path: Path):
    """should_skip_step must return False and reset status if recorded output file is deleted."""
    state = PipelineState(output_dir=str(tmp_path), enable_checksum=False)
    state._loaded_from_file = True

    out_file = tmp_path / "stage1.out"
    out_file.write_text("hello world\n")

    state.start_step("stage1")
    state.complete_step("stage1", output_files=[str(out_file)])
    assert state.state["steps"]["stage1"].status == "completed"

    # With file existing, should_skip_step returns True
    assert state.should_skip_step("stage1") is True

    # Delete the output file - should_skip_step must now return False and reset status to pending
    out_file.unlink()
    assert state.should_skip_step("stage1") is False
    assert state.state["steps"]["stage1"].status == "pending"


@pytest.mark.unit
def test_should_skip_step_validates_output_files_checksum(tmp_path: Path):
    """should_skip_step must return False if checksum mismatch occurs with enable_checksum=True."""
    state = PipelineState(output_dir=str(tmp_path), enable_checksum=True)
    state._loaded_from_file = True

    out_file = tmp_path / "stage1.out"
    out_file.write_text("original content")

    state.start_step("stage1")
    state.complete_step("stage1", output_files=[str(out_file)])
    assert state.state["steps"]["stage1"].status == "completed"

    # Modify file content with same or different size
    out_file.write_text("tampered content")
    assert state.should_skip_step("stage1") is False
    assert state.state["steps"]["stage1"].status == "pending"


@pytest.mark.unit
def test_pipeline_runner_should_skip_stage_method(tmp_path: Path):
    """PipelineRunner must provide _should_skip_stage method."""
    runner = PipelineRunner(enable_checkpoints=True, max_workers=1)
    context = PipelineContext(args=Mock(), config={"output_dir": str(tmp_path)}, workspace=Mock())
    stage = DummyStage("test_stage")

    assert hasattr(runner, "_should_skip_stage")
    assert callable(runner._should_skip_stage)
    # When not completed, should return False
    assert runner._should_skip_stage(stage, context) is False


@pytest.mark.unit
def test_upstream_rerun_invalidates_downstream_stages(tmp_path: Path):
    """When an upstream stage re-runs, downstream stages must be invalidated and re-executed."""
    state = PipelineState(output_dir=str(tmp_path), enable_checksum=False)
    state.initialize({"output_dir": str(tmp_path), "resume": True}, pipeline_version="1.0")

    context = PipelineContext(
        args=Mock(),
        config={"output_dir": str(tmp_path), "resume": True, "pipeline_version": "1.0"},
        workspace=Mock(),
    )
    context.checkpoint_state = state

    s1 = DummyStage("s1")
    s2 = DummyStage("s2", dependencies=["s1"])
    s3 = DummyStage("s3", dependencies=["s2"])
    stages = [s1, s2, s3]

    runner = PipelineRunner(enable_checkpoints=True, max_workers=1)

    # First run: all 3 stages run
    runner.run(stages, context)
    assert s1.run_count == 1
    assert s2.run_count == 1
    assert s3.run_count == 1

    # Second run without changes: all 3 stages are skipped
    c2 = PipelineContext(
        args=Mock(),
        config={"output_dir": str(tmp_path), "resume": True, "pipeline_version": "1.0"},
        workspace=Mock(),
    )
    s1_b = DummyStage("s1")
    s2_b = DummyStage("s2", dependencies=["s1"])
    s3_b = DummyStage("s3", dependencies=["s2"])

    state2 = PipelineState(output_dir=str(tmp_path), enable_checksum=False)
    state2.load()
    c2.checkpoint_state = state2

    runner.run([s1_b, s2_b, s3_b], c2)
    assert s1_b.run_count == 0
    assert s2_b.run_count == 0
    assert s3_b.run_count == 0

    # Now delete s1's output file!
    (tmp_path / "s1.out").unlink()

    c3 = PipelineContext(
        args=Mock(),
        config={"output_dir": str(tmp_path), "resume": True, "pipeline_version": "1.0"},
        workspace=Mock(),
    )
    s1_c = DummyStage("s1")
    s2_c = DummyStage("s2", dependencies=["s1"])
    s3_c = DummyStage("s3", dependencies=["s2"])

    state3 = PipelineState(output_dir=str(tmp_path), enable_checksum=False)
    state3.load()
    c3.checkpoint_state = state3

    runner.run([s1_c, s2_c, s3_c], c3)
    # Because s1 was re-run, downstream stages s2 and s3 must also re-run!
    assert s1_c.run_count == 1
    assert s2_c.run_count == 1
    assert s3_c.run_count == 1


@pytest.mark.unit
def test_stage_enable_checkpoint_flag_disabled(tmp_path: Path):
    """When a stage has enable_checkpoint=False, it is never skipped by checkpointing."""
    state = PipelineState(output_dir=str(tmp_path))
    state._loaded_from_file = True
    state.start_step("no_chk")
    state.complete_step("no_chk")

    stage = DummyStage("no_chk", enable_chk=False)
    context = PipelineContext(args=Mock(), config={"output_dir": str(tmp_path)}, workspace=Mock())
    context.checkpoint_state = state

    runner = PipelineRunner(enable_checkpoints=True, max_workers=1)
    assert runner._should_skip_stage(stage, context) is False

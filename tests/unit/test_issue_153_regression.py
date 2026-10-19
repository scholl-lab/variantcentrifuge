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
    state.initialize(
        {"output_dir": str(tmp_path), "enable_checkpoint": True}, pipeline_version="1.0"
    )

    context = PipelineContext(
        args=Mock(),
        config={"output_dir": str(tmp_path), "enable_checkpoint": True, "pipeline_version": "1.0"},
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


@pytest.mark.unit
@pytest.mark.parametrize(
    "key,val1,val2",
    [
        ("filter_exception_file", "hotspots_a.tsv", "hotspots_b.tsv"),
        ("transcript_list", ["NM_001"], ["NM_002"]),
        ("transcript_file", "tx_a.txt", "tx_b.txt"),
        ("regions_bed", "reg_a.bed", "reg_b.bed"),
        ("interval_expand", 50, 100),
        ("add_chr", True, False),
        ("tumor_sample_index", "0", "1"),
        ("normal_sample_index", "1", "2"),
    ],
)
def test_changing_named_configuration_keys_invalidates_compatibility(
    tmp_path: Path, key, val1, val2
):
    """Changing each of the 8 named configuration keys invalidates checkpoint compatibility."""
    from variantcentrifuge.checkpoint_helpers import hash_configuration

    cfg_base = {"output_dir": str(tmp_path), "vcf_file": "input.vcf.gz", "gene_name": "GENE1"}
    cfg1 = {**cfg_base, key: val1}
    cfg2 = {**cfg_base, key: val2}

    h1 = hash_configuration(cfg1)
    h2 = hash_configuration(cfg2)
    assert h1 != h2, f"Hash did not change when modifying {key}"

    state = PipelineState(output_dir=str(tmp_path))
    state.initialize(cfg1, pipeline_version="1.0")
    state._loaded_from_file = True

    assert state.can_resume(cfg1, pipeline_version="1.0") is True
    assert state.can_resume(cfg2, pipeline_version="1.0") is False


@pytest.mark.unit
def test_legacy_fingerprint_fallback_rejected_when_omitted_parameter_changes(tmp_path: Path):
    """Older unversioned fingerprint cannot regain compatibility through legacy fallback
    when omitted parameter is used or removed, and unversioned v1 hashes are rejected.
    """
    from variantcentrifuge.checkpoint_helpers import hash_legacy_configuration

    # Case 1: unversioned checkpoint without transcript_list cannot resume with new or old config
    cfg_old = {"output_dir": str(tmp_path), "vcf_file": "input.vcf.gz"}
    legacy_hash = hash_legacy_configuration(cfg_old)

    state = PipelineState(output_dir=str(tmp_path))
    state.initialize(cfg_old, pipeline_version="1.0")
    state.state["configuration_hash"] = legacy_hash
    state._loaded_from_file = True

    assert state.can_resume(cfg_old, pipeline_version="1.0") is False
    cfg_new = {**cfg_old, "transcript_list": ["NM_NEW"]}
    assert state.can_resume(cfg_new, pipeline_version="1.0") is False

    # Case 2: unversioned checkpoint with transcript_list cannot resume when transcript_list removed
    cfg_old_with_tx = {
        "output_dir": str(tmp_path),
        "vcf_file": "input.vcf.gz",
        "transcript_list": ["NM_OLD"],
    }
    legacy_hash_tx = hash_legacy_configuration(cfg_old_with_tx)
    state2 = PipelineState(output_dir=str(tmp_path))
    state2.initialize(cfg_old_with_tx, pipeline_version="1.0")
    state2.state["configuration_hash"] = legacy_hash_tx
    state2._loaded_from_file = True

    cfg_new_without_tx = {"output_dir": str(tmp_path), "vcf_file": "input.vcf.gz"}
    assert state2.can_resume(cfg_new_without_tx, pipeline_version="1.0") is False


@pytest.mark.unit
def test_incompatible_traditional_resume_stops_immediately(tmp_path: Path):
    """Incompatible traditional resume must stop/abort immediately without stage skipping."""
    state = PipelineState(output_dir=str(tmp_path))
    state.initialize(
        {"output_dir": str(tmp_path), "transcript_list": ["NM_OLD"]},
        pipeline_version="1.0",
    )
    state.start_step("s1")
    out_file = tmp_path / "s1.out"
    out_file.write_text("old data")
    state.complete_step("s1", output_files=[str(out_file)])
    state.save()

    new_cfg = {
        "output_dir": str(tmp_path),
        "transcript_list": ["NM_NEW"],
        "resume": True,
        "pipeline_version": "1.0",
    }
    context = PipelineContext(args=Mock(), config=new_cfg, workspace=Mock())
    state2 = PipelineState(output_dir=str(tmp_path))
    state2.load()
    context.checkpoint_state = state2

    s1 = DummyStage("s1")
    runner = PipelineRunner(enable_checkpoints=True, max_workers=1)

    with pytest.raises(ValueError, match="Cannot resume: configuration or version mismatch"):
        runner.run([s1], context)

    assert s1.run_count == 0


@pytest.mark.unit
def test_registered_stages_resolve_with_correct_categories():
    """Verify PCAComputationStage, DataSortingStage, TranscriptFilterStage,
    and ClinVarPM5Stage in registry.
    """
    from variantcentrifuge.stages.stage_registry import get_registry, initialize_registry

    initialize_registry()
    registry = get_registry()

    pca_info = registry.get_stage_info("pca_computation")
    assert pca_info is not None
    assert pca_info.category == "processing"

    sort_info = registry.get_stage_info("data_sorting")
    assert sort_info is not None
    assert sort_info.category == "processing"

    tx_info = registry.get_stage_info("transcript_filter")
    assert tx_info is not None
    assert tx_info.category == "processing"

    pm5_info = registry.get_stage_info("clinvar_pm5")
    assert pm5_info is not None
    assert pm5_info.category == "analysis"


@pytest.mark.unit
def test_completing_tsv_output_alone_does_not_bypass_artifact_validation(tmp_path: Path):
    """Completing tsv_output alone does not bypass validation of required
    artifacts or subsequent reports.
    """
    state = PipelineState(output_dir=str(tmp_path))
    cfg = {"output_dir": str(tmp_path), "vcf_file": "in.vcf.gz"}
    state.initialize(cfg, pipeline_version="1.0")
    state._loaded_from_file = True

    # Intermediate stage with recorded output file that is deleted
    missing_intermediate = tmp_path / "intermediate.tsv"
    missing_intermediate.write_text("intermediate data")
    state.start_step("field_extraction")
    state.complete_step("field_extraction", output_files=[str(missing_intermediate)])
    missing_intermediate.unlink()

    # tsv_output completed with existing file
    tsv_out = tmp_path / "final.tsv"
    tsv_out.write_text("complete final tsv")
    state.start_step("tsv_output")
    state.complete_step("tsv_output", output_files=[str(tsv_out)])

    # can_resume must fail because intermediate artifact is missing, even though tsv_output exists
    assert state.can_resume(cfg, pipeline_version="1.0") is False

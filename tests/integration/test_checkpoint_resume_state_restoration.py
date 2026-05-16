"""Regression tests for checkpoint resume state restoration."""

from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.checkpoint import PipelineState
from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.runner import PipelineRunner
from variantcentrifuge.pipeline_core.stage import Stage
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.analysis_stages import ChunkedAnalysisStage, DataFrameLoadingStage


class VolatileAnnotationStage(Stage):
    @property
    def name(self) -> str:
        return "custom_annotation"

    @property
    def dependencies(self) -> set[str]:
        return {"dataframe_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        return {"chunked_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        assert context.current_dataframe is not None
        context.current_dataframe["RESTORED_MARKER"] = "yes"
        return context


class AssociationProbeStage(Stage):
    def __init__(self, output_name: str = "association.tsv"):
        super().__init__()
        self.saw_dataframe = False
        self.output_name = output_name

    @property
    def name(self) -> str:
        return "association_analysis"

    @property
    def dependencies(self) -> set[str]:
        return {"dataframe_loading", "custom_annotation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        assert context.current_dataframe is not None
        assert "RESTORED_MARKER" in context.current_dataframe.columns
        self.saw_dataframe = True
        assoc_output = context.workspace.output_dir / self.output_name
        pd.DataFrame(columns=["gene", "primary_pvalue"]).to_csv(assoc_output, sep="\t", index=False)
        context.config["association_output"] = str(assoc_output)
        return context


class FinalTSVStage(Stage):
    @property
    def name(self) -> str:
        return "tsv_output"

    @property
    def dependencies(self) -> set[str]:
        return {"association_analysis"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        assert context.current_dataframe is not None
        output = context.workspace.output_dir / "variants.tsv"
        context.current_dataframe.to_csv(output, sep="\t", index=False)
        context.final_output_path = output
        return context


def _context_with_state(tmp_path: Path, chunked_text: str | None) -> PipelineContext:
    workspace = Workspace(tmp_path, "resume_case")
    workspace.intermediate_dir.mkdir(parents=True, exist_ok=True)
    extracted = workspace.intermediate_dir / "resume_case.extracted.tsv"
    extracted.write_text("GENE\tCHROM\nBRCA1\tchr1\n")
    if chunked_text is not None:
        (workspace.intermediate_dir / "chunked_analysis_results.tsv").write_text(chunked_text)

    config = {
        "resume": True,
        "pipeline_version": "test-version",
        "output_file": "variants.tsv",
        "perform_association": True,
        "gzip_intermediates": False,
        "force_chunked_processing": True,
    }
    state = PipelineState(str(tmp_path), enable_checksum=False)
    state.initialize(config, "test-version")
    for name in [
        "dataframe_loading",
        "chunked_analysis",
        "custom_annotation",
        "association_analysis",
        "tsv_output",
    ]:
        state.start_step(name)
        state.complete_step(name, [], [])

    resume_state = PipelineState(str(tmp_path), enable_checksum=False)
    assert resume_state.load()
    context = PipelineContext(args=Namespace(), config=config, workspace=workspace)
    context.data = extracted
    context.extracted_tsv = extracted
    context.checkpoint_state = resume_state
    return context


def test_plain_resume_restores_chunked_dataframe_before_association(tmp_path):
    """Issue #101 shape: chunked_analysis restore happens before association resumes."""
    context = _context_with_state(tmp_path, "GENE\tCHROM\nBRCA1\tchr1\n")
    association = AssociationProbeStage()
    stages = [
        DataFrameLoadingStage(),
        ChunkedAnalysisStage(),
        VolatileAnnotationStage(),
        association,
        FinalTSVStage(),
    ]

    result = PipelineRunner(enable_checkpoints=True, max_workers=1).run(stages, context)

    assert association.saw_dataframe is True
    assert result.current_dataframe is not None
    assert result.final_output_path.exists()
    assert Path(result.config["association_output"]).exists()
    assert "RESTORED_MARKER" in result.final_output_path.read_text()


def test_plain_resume_aborts_when_chunked_artifact_missing(tmp_path):
    """Missing chunked analysis artifact must abort resume before downstream success."""
    context = _context_with_state(tmp_path, None)
    stages = [
        DataFrameLoadingStage(),
        ChunkedAnalysisStage(),
        VolatileAnnotationStage(),
        AssociationProbeStage(),
        FinalTSVStage(),
    ]

    with pytest.raises(RuntimeError, match="Cannot restore chunked_analysis"):
        PipelineRunner(enable_checkpoints=True, max_workers=1).run(stages, context)


def test_plain_resume_restores_header_only_chunked_dataframe(tmp_path):
    """Header-only chunked output is a valid empty resumed DataFrame."""
    context = _context_with_state(tmp_path, "GENE\tCHROM\n")
    stages = [
        DataFrameLoadingStage(),
        ChunkedAnalysisStage(),
        VolatileAnnotationStage(),
        AssociationProbeStage(),
        FinalTSVStage(),
    ]

    result = PipelineRunner(enable_checkpoints=True, max_workers=1).run(stages, context)

    assert result.current_dataframe is not None
    assert result.current_dataframe.empty
    assert result.final_output_path.read_text().startswith("GENE\tCHROM\tRESTORED_MARKER\n")

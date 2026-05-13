from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.pipeline_core.context import PipelineContext
from variantcentrifuge.pipeline_core.workspace import Workspace
from variantcentrifuge.stages.processing_stages import TranscriptFilterStage


@pytest.fixture
def context(tmp_path: Path) -> PipelineContext:
    workspace = Workspace(output_dir=tmp_path, base_name="sample")
    return PipelineContext(
        args=None,
        config={},
        workspace=workspace,
        data=None,
    )


def test_transcript_filter_skips_without_requested_transcripts(
    context: PipelineContext,
) -> None:
    input_vcf = Path("input.vcf.gz")
    context.extracted_vcf = input_vcf
    context.data = input_vcf
    context.config = {}

    result = TranscriptFilterStage()._process(context)

    assert result.data == input_vcf
    assert not hasattr(result, "transcript_filtered_vcf")


@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
def test_transcript_filter_uses_extracted_vcf_and_sets_context(
    mock_filter,
    context: PipelineContext,
    tmp_path: Path,
) -> None:
    input_vcf = tmp_path / "input.vcf.gz"
    input_vcf.write_text("dummy", encoding="utf-8")
    context.extracted_vcf = input_vcf
    context.data = input_vcf
    context.config = {"transcript_list": "NM_000059.4,NM_007294.4"}

    result = TranscriptFilterStage()._process(context)

    mock_filter.assert_called_once()
    assert mock_filter.call_args.args[0] == str(input_vcf)
    assert mock_filter.call_args.args[2] == {"NM_000059.4", "NM_007294.4"}
    assert result.transcript_filtered_vcf.name == "sample.transcript_filtered.vcf.gz"
    assert result.data == result.transcript_filtered_vcf

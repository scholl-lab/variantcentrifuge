from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.stages.processing_stages import ParallelCompleteProcessingStage


@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_applies_transcript_filter_before_field_extraction(
    mock_extract_variants,
    mock_filter_transcripts,
    mock_extract_fields,
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    stage = ParallelCompleteProcessingStage()
    result, _ = stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "fields_to_extract": "CHROM POS ANN[0].FEATUREID GEN[*].GT",
            "transcript_ids": {"NM_keep.1"},
            "vcf_samples": ["S1"],
        },
    )

    transcript_vcf = tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    mock_filter_transcripts.assert_called_once()
    assert mock_filter_transcripts.call_args.args[2] == {"NM_keep.1"}
    assert mock_filter_transcripts.call_args.kwargs["index_output"] is False
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(transcript_vcf)
    assert result.name == "sample.chunk_0.extracted.tsv.gz"


@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.load_transcript_ids")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_uses_preloaded_transcript_ids(
    mock_extract_variants,
    mock_filter_transcripts,
    mock_load_transcript_ids,
    mock_extract_fields,
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    stage = ParallelCompleteProcessingStage()
    stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "fields_to_extract": "CHROM POS ANN[0].FEATUREID GEN[*].GT",
            "transcript_ids": {"NM_parent_loaded.1"},
            "transcript_file": "transcripts.txt",
            "vcf_samples": ["S1"],
        },
    )

    mock_load_transcript_ids.assert_not_called()
    mock_filter_transcripts.assert_called_once()
    assert mock_filter_transcripts.call_args.args[2] == {"NM_parent_loaded.1"}


@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
@patch("variantcentrifuge.stages.processing_stages.split_snpeff_annotations")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_splits_before_snpsift_when_requested(
    mock_extract_variants,
    mock_split,
    mock_apply_filter,
    mock_filter_transcripts,
    mock_extract_fields,
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    stage = ParallelCompleteProcessingStage()
    stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "snpeff_splitting_mode": "before_filters",
            "filters": "(ANN[ANY].IMPACT has 'HIGH')",
            "fields_to_extract": "CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
            "transcript_ids": {"NM_gene_b.1"},
            "vcf_samples": ["S1"],
        },
    )

    chunk_vcf = tmp_path / "sample.chunk_0.variants.vcf.gz"
    split_vcf = tmp_path / "sample.chunk_0.split_annotations.vcf.gz"
    transcript_vcf = tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    filtered_vcf = tmp_path / "sample.chunk_0.filtered.vcf.gz"

    mock_split.assert_called_once_with(str(chunk_vcf), str(split_vcf))
    assert mock_filter_transcripts.call_args.args[0] == str(split_vcf)
    assert mock_apply_filter.call_args.args[0] == str(transcript_vcf)
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(filtered_vcf)


@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
@patch("variantcentrifuge.stages.processing_stages.split_snpeff_annotations")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_splits_after_snpsift_when_requested(
    mock_extract_variants,
    mock_split,
    mock_apply_filter,
    mock_filter_transcripts,
    mock_extract_fields,
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    stage = ParallelCompleteProcessingStage()
    stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "snpeff_splitting_mode": "after_filters",
            "filters": "(ANN[ANY].IMPACT has 'HIGH')",
            "fields_to_extract": "CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
            "transcript_ids": {"NM_gene_b.1"},
            "vcf_samples": ["S1"],
        },
    )

    chunk_vcf = tmp_path / "sample.chunk_0.variants.vcf.gz"
    transcript_vcf = tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    filtered_vcf = tmp_path / "sample.chunk_0.filtered.vcf.gz"
    split_vcf = tmp_path / "sample.chunk_0.split_annotations.vcf.gz"

    assert mock_filter_transcripts.call_args.args[0] == str(chunk_vcf)
    assert mock_apply_filter.call_args.args[0] == str(transcript_vcf)
    mock_split.assert_called_once_with(str(filtered_vcf), str(split_vcf))
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(split_vcf)


@pytest.mark.parametrize(
    "extra_config",
    [
        {"late_filtering": True, "filters": "(ANN[ANY].IMPACT has 'HIGH')"},
        {},
    ],
)
@patch("variantcentrifuge.stages.processing_stages.extract_fields_bcftools")
@patch("variantcentrifuge.stages.processing_stages.filter_vcf_to_transcripts")
@patch("variantcentrifuge.stages.processing_stages.apply_snpsift_filter")
@patch("variantcentrifuge.stages.processing_stages.split_snpeff_annotations")
@patch("variantcentrifuge.stages.processing_stages.extract_variants")
def test_parallel_chunk_before_filter_fallbacks_preserve_split_working_vcf(
    mock_extract_variants,
    mock_split,
    mock_apply_filter,
    mock_filter_transcripts,
    mock_extract_fields,
    extra_config: dict[str, object],
    tmp_path: Path,
) -> None:
    chunk_bed = tmp_path / "chunk.bed"
    chunk_bed.write_text("1\t0\t100\n", encoding="utf-8")

    config = {
        "threads_per_chunk": 1,
        "snpeff_splitting_mode": "before_filters",
        "fields_to_extract": "CHROM POS ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
        "transcript_ids": {"NM_gene_b.1"},
        "vcf_samples": ["S1"],
    }
    config.update(extra_config)

    stage = ParallelCompleteProcessingStage()
    stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config=config,
    )

    chunk_vcf = tmp_path / "sample.chunk_0.variants.vcf.gz"
    split_vcf = tmp_path / "sample.chunk_0.split_annotations.vcf.gz"

    mock_split.assert_called_once_with(str(chunk_vcf), str(split_vcf))
    mock_apply_filter.assert_not_called()
    assert mock_filter_transcripts.call_args.args[0] == str(split_vcf)
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(
        tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    )


@pytest.mark.parametrize("split_mode", [None, "before_filters", "after_filters"])
def test_sequential_pipeline_stage_order_with_transcripts_and_split_modes(
    split_mode: str | None,
) -> None:
    """Verify build_pipeline_stages + runner.dry_run order for all split modes with transcripts."""
    from variantcentrifuge.cli import create_parser
    from variantcentrifuge.pipeline import build_pipeline_stages
    from variantcentrifuge.pipeline_core import PipelineRunner

    parser = create_parser()
    cli_argv = [
        "-v",
        "input.vcf.gz",
        "-g",
        "TP53",
        "-f",
        "(ANN[ANY].IMPACT has 'HIGH')",
        "--transcript-list",
        "NM_000546.6",
        "--threads",
        "1",
    ]
    if split_mode is not None:
        cli_argv.extend(["--split-snpeff-lines", split_mode])
    args = parser.parse_args(cli_argv)
    args.threads = 1
    args.config = {
        "snpeff_splitting_mode": split_mode,
        "transcript_list": "NM_000546.6",
        "filters": "(ANN[ANY].IMPACT has 'HIGH')",
    }

    stages = build_pipeline_stages(args)
    plan = PipelineRunner().dry_run(stages)
    level_of = {name: lvl for lvl, names in enumerate(plan) for name in names}

    assert level_of["variant_extraction"] < level_of["transcript_filtering"]
    assert level_of["transcript_filtering"] < level_of["snpsift_filtering"]
    assert level_of["snpsift_filtering"] < level_of["field_extraction"]
    if split_mode == "before_filters":
        assert level_of["multiallelic_split"] < level_of["transcript_filtering"]
        assert level_of["multiallelic_split"] < level_of["snpsift_filtering"]
    elif split_mode == "after_filters":
        assert level_of["snpsift_filtering"] < level_of["multiallelic_split"]
        assert level_of["multiallelic_split"] < level_of["field_extraction"]

from pathlib import Path
from unittest.mock import patch

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
    result = stage._process_single_chunk(
        chunk_index=0,
        chunk_bed=chunk_bed,
        vcf_file="input.vcf.gz",
        base_name="sample",
        intermediate_dir=tmp_path,
        config={
            "threads_per_chunk": 1,
            "fields_to_extract": "CHROM POS ANN[0].FEATUREID GEN[*].GT",
            "transcript_list": "NM_keep.1",
            "vcf_samples": ["S1"],
        },
    )

    transcript_vcf = tmp_path / "sample.chunk_0.transcript_filtered.vcf.gz"
    mock_filter_transcripts.assert_called_once()
    assert mock_filter_transcripts.call_args.args[2] == {"NM_keep.1"}
    assert mock_extract_fields.call_args.kwargs["variant_file"] == str(transcript_vcf)
    assert result.name == "sample.chunk_0.extracted.tsv.gz"

import gzip
from pathlib import Path
from unittest.mock import patch

import pytest

from variantcentrifuge.transcript_filter import (
    filter_ann_value_by_transcripts,
    filter_vcf_line_by_transcripts,
    filter_vcf_to_transcripts,
    load_transcript_ids,
)


def test_load_transcript_ids_from_list_and_file(tmp_path: Path) -> None:
    transcript_file = tmp_path / "transcripts.txt"
    transcript_file.write_text(
        "# comment\nNM_000059.4\n\nENST00000357654.9\n",
        encoding="utf-8",
    )

    result = load_transcript_ids("NM_007294.4, NM_000546.6", str(transcript_file))

    assert result == {
        "NM_007294.4",
        "NM_000546.6",
        "NM_000059.4",
        "ENST00000357654.9",
    }


def test_load_transcript_ids_rejects_empty_file(tmp_path: Path) -> None:
    transcript_file = tmp_path / "empty.txt"
    transcript_file.write_text("# only comments\n\n", encoding="utf-8")

    with pytest.raises(ValueError, match="No transcript IDs"):
        load_transcript_ids(None, str(transcript_file))


def test_filter_ann_keeps_matching_feature_id_and_reorders_to_single_ann() -> None:
    ann = (
        "A|synonymous_variant|LOW|BRCA2|BRCA2|transcript|NM_other.1|protein_coding||||||||,"
        "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_000059.4|protein_coding|"
        "10/27|c.1234G>A|p.Arg412His|||||,"
        "A|intron_variant|MODIFIER|BRCA2|BRCA2|transcript|ENST00000380152.8|"
        "protein_coding||||||||"
    )

    result = filter_ann_value_by_transcripts(ann, {"NM_000059.4"})

    assert result == (
        "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_000059.4|protein_coding|"
        "10/27|c.1234G>A|p.Arg412His|||||"
    )


def test_filter_ann_keeps_all_matching_feature_ids_in_original_order() -> None:
    ann = (
        "A|missense_variant|MODERATE|GENE|GENE|transcript|NM_1|protein_coding||||||||,"
        "A|splice_region_variant|LOW|GENE|GENE|transcript|NM_2|protein_coding||||||||,"
        "A|intron_variant|MODIFIER|GENE|GENE|transcript|NM_3|protein_coding||||||||"
    )

    result = filter_ann_value_by_transcripts(ann, {"NM_1", "NM_3"})

    assert result == (
        "A|missense_variant|MODERATE|GENE|GENE|transcript|NM_1|protein_coding||||||||,"
        "A|intron_variant|MODIFIER|GENE|GENE|transcript|NM_3|protein_coding||||||||"
    )


def test_filter_ann_returns_none_when_no_transcript_matches() -> None:
    ann = "A|missense_variant|MODERATE|BRCA2|BRCA2|transcript|NM_other.1|protein_coding||||||||"

    assert filter_ann_value_by_transcripts(ann, {"NM_000059.4"}) is None


def test_filter_vcf_line_rewrites_ann_and_preserves_other_info() -> None:
    line = (
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "AC=1;ANN=G|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|"
        "protein_coding||||||||,"
        "G|missense_variant|MODERATE|GENE|GENE|transcript|NM_keep.1|"
        "protein_coding||||||||;DP=20"
        "\tGT\t0/1\n"
    )

    result = filter_vcf_line_by_transcripts(line, {"NM_keep.1"})

    assert result is not None
    assert "AC=1;" in result
    assert "DP=20" in result
    assert "NM_keep.1" in result
    assert "NM_other.1" not in result


def test_filter_vcf_line_drops_no_match_record() -> None:
    line = (
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "ANN=G|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|"
        "protein_coding||||||||"
        "\tGT\t0/1\n"
    )

    assert filter_vcf_line_by_transcripts(line, {"NM_keep.1"}) is None


def test_filter_vcf_to_transcripts_writes_only_matching_records(tmp_path: Path) -> None:
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf.gz"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "ANN=G|missense_variant|MODERATE|GENE|GENE|transcript|NM_keep.1|"
        "protein_coding||||||||\tGT\t0/1\n"
        "1\t200\t.\tC\tT\t.\tPASS\t"
        "ANN=T|synonymous_variant|LOW|GENE|GENE|transcript|NM_other.1|"
        "protein_coding||||||||\tGT\t0/1\n",
        encoding="utf-8",
    )

    count = filter_vcf_to_transcripts(str(input_vcf), str(output_vcf), {"NM_keep.1"})

    assert count == 1
    with gzip.open(output_vcf, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "#CHROM" in text
    assert "NM_keep.1" in text
    assert "NM_other.1" not in text


@patch("variantcentrifuge.transcript_filter.subprocess.run")
def test_filter_vcf_to_transcripts_uses_bcftools_bgzip_when_indexing(
    mock_run,
    tmp_path: Path,
) -> None:
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf.gz"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "1\t100\t.\tA\tG\t.\tPASS\t"
        "ANN=G|missense_variant|MODERATE|GENE|GENE|transcript|NM_keep.1|"
        "protein_coding||||||||\tGT\t0/1\n",
        encoding="utf-8",
    )

    count = filter_vcf_to_transcripts(
        str(input_vcf),
        str(output_vcf),
        {"NM_keep.1"},
        index_output=True,
    )

    assert count == 1
    assert mock_run.call_count == 2
    view_call = mock_run.call_args_list[0].args[0]
    index_call = mock_run.call_args_list[1].args[0]
    assert view_call[:4] == ["bcftools", "view", "-Oz", "-o"]
    assert view_call[4] == str(output_vcf)
    assert index_call == ["bcftools", "index", "-f", str(output_vcf)]

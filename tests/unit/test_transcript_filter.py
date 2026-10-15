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


def test_filter_vcf_to_transcripts_matches_version_skew_emits_each_kept_and_logs_summary(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        "chr3\t53820000\t.\tT\tC\t.\tPASS\t"
        "ANN=C|stop_gained|HIGH|OTHER|OTHER|transcript|NM_999999.1|protein_coding||||||||,"
        "C|missense_variant|MODERATE|CACNA1D|CACNA1D|transcript|NM_000720.3|protein_coding|"
        "8/49|c.1207G>C|p.Gly403Arg|||||,"
        "C|missense_variant|MODERATE|CACNA1D|CACNA1D|transcript|NM_001128840.3|protein_coding|"
        "8/48|c.1207G>C|p.Gly403Arg|||||\tGT\t0/1\n",
        encoding="utf-8",
    )

    with caplog.at_level("INFO", logger="variantcentrifuge.transcript_filter"):
        count = filter_vcf_to_transcripts(
            str(input_vcf),
            str(output_vcf),
            {"NM_000720.4", "NM_001128840", "NM_unmatched.2"},
        )

    assert count == 2
    lines = [
        line
        for line in output_vcf.read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    assert len(lines) == 2
    assert "NM_000720.3" in lines[0] and "NM_001128840.3" not in lines[0]
    assert "NM_001128840.3" in lines[1] and "NM_000720.3" not in lines[1]
    assert "NM_999999.1" not in output_vcf.read_text(encoding="utf-8")
    assert "matched 2/3" in caplog.text
    assert "NM_unmatched.2" in caplog.text
    assert "NM_000720.4->NM_000720.3" in caplog.text


def test_filter_exception_file_rescues_matching_hotspot_only(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Hotspot with FILTER=normal_artifact is rescued, while strand_bias or wrong ALT is not."""
    from variantcentrifuge.filters import apply_filter_exceptions_to_vcf, load_filter_exceptions

    exc_file = tmp_path / "exceptions.tsv"
    exc_file.write_text(
        "CHROM\tPOS\tREF\tALT\tALLOWED_FILTERS\n"
        "chr11\t128839188\tG\tA\tnormal_artifact,panel_of_normals\n",
        encoding="utf-8",
    )
    in_vcf = tmp_path / "in.vcf"
    out_vcf = tmp_path / "out.vcf"
    in_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n"
        "chr11\t128839188\t.\tG\tA\t.\tnormal_artifact\tDP=100\tGT\t0/0\t0/1\n"
        "chr11\t128839188\t.\tG\tA\t.\tnormal_artifact;strand_bias\tDP=100\tGT\t0/0\t0/1\n"
        "chr11\t128839188\t.\tG\tT\t.\tnormal_artifact\tDP=100\tGT\t0/0\t0/1\n",
        encoding="utf-8",
    )

    exceptions = load_filter_exceptions(str(exc_file))
    with caplog.at_level("INFO", logger="variantcentrifuge"):
        rescued = apply_filter_exceptions_to_vcf(str(in_vcf), str(out_vcf), exceptions)

    assert rescued == 1
    assert (
        "Rescued hotspot variant chr11:128839188 G>A (RAW_FILTER=normal_artifact -> PASS)"
        in caplog.text
    )
    records = [
        line.split("\t")
        for line in out_vcf.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert records[0][6] == "PASS"
    assert "RAW_FILTER=normal_artifact" in records[0][7]
    assert records[1][6] == "normal_artifact;strand_bias"
    assert "RAW_FILTER=" not in records[1][7]
    assert records[2][6] == "normal_artifact"
    assert "RAW_FILTER=" not in records[2][7]

    # Verify filter_vcf_to_transcripts emits ##INFO=<ID=RAW_FILTER,...> header
    tx_in = tmp_path / "tx_in.vcf"
    tx_out = tmp_path / "tx_out.vcf"
    tx_in.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tN\tT\n"
        "chr11\t128839188\t.\tG\tA\t.\tnormal_artifact\t"
        "ANN=A|missense_variant|MODERATE|KCNJ5|KCNJ5|transcript|NM_000890.5|protein_coding||||||||"
        "\tGT\t0/0\t0/1\n",
        encoding="utf-8",
    )
    filter_vcf_to_transcripts(
        str(tx_in),
        str(tx_out),
        {"NM_000890.5"},
        filter_exception_file=str(exc_file),
    )
    tx_out_text = tx_out.read_text(encoding="utf-8")
    assert "##INFO=<ID=RAW_FILTER,Number=1,Type=String" in tx_out_text
    assert "RAW_FILTER=normal_artifact" in tx_out_text


def test_exact_version_match_preferred_over_base_id_fallback() -> None:
    """When both NM_000720.1 and NM_000720.4 appear in ANN, requesting .4 keeps only .4."""
    ann = (
        "A|synonymous_variant|LOW|CACNA1D|CACNA1D|transcript|NM_000720.1|protein_coding||||||||,"
        "A|missense_variant|MODERATE|CACNA1D|CACNA1D|transcript|NM_000720.4|protein_coding||||||||"
    )
    result = filter_ann_value_by_transcripts(ann, {"NM_000720.4"})
    assert result == (
        "A|missense_variant|MODERATE|CACNA1D|CACNA1D|transcript|NM_000720.4|protein_coding||||||||"
    )

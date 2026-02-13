"""Unit tests for variantcentrifuge.vcf_header_parser."""

from __future__ import annotations

import json
from unittest.mock import patch

import pytest

from variantcentrifuge.vcf_header_parser import (
    AnnotationField,
    format_fields_json,
    format_fields_table,
    parse_vcf_header,
    validate_requested_fields,
)

# ---------------------------------------------------------------------------
# Fixture: realistic VCF header text
# ---------------------------------------------------------------------------
SAMPLE_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation'">
##INFO=<ID=dbNSFP_REVEL_score,Number=.,Type=String,Description="REVEL score">
##INFO=<ID=dbNSFP_CADD_phred,Number=.,Type=String,Description="CADD phred-like score">
##INFO=<ID=dbNSFP_GERP++_RS,Number=.,Type=String,Description="GERP++ RS score">
##INFO=<ID=dbNSFP_hg19_chr,Number=.,Type=String,Description="hg19 chromosome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


@pytest.fixture()
def sample_fields() -> list[AnnotationField]:
    """Parse the sample header and return fields."""
    with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=SAMPLE_HEADER):
        return parse_vcf_header("dummy.vcf")


# ---------------------------------------------------------------------------
# parse_vcf_header
# ---------------------------------------------------------------------------
class TestParseVcfHeader:
    def test_parses_info_fields(self, sample_fields: list[AnnotationField]) -> None:
        info = [f for f in sample_fields if f.source == "INFO"]
        assert len(info) == 7
        ids = {f.id for f in info}
        assert "AC" in ids
        assert "dbNSFP_REVEL_score" in ids
        assert "dbNSFP_GERP++_RS" in ids

    def test_parses_format_fields(self, sample_fields: list[AnnotationField]) -> None:
        fmt = [f for f in sample_fields if f.source == "FORMAT"]
        assert len(fmt) == 4
        ids = {f.id for f in fmt}
        assert "GT" in ids
        assert "DP" in ids

    def test_sorted_by_source_and_id(self, sample_fields: list[AnnotationField]) -> None:
        info = [f for f in sample_fields if f.source == "INFO"]
        format_fields = [f for f in sample_fields if f.source == "FORMAT"]
        # Alphabetically: FORMAT < INFO, so FORMAT fields come first
        if info and format_fields:
            assert max(sample_fields.index(f) for f in format_fields) < min(
                sample_fields.index(f) for f in info
            )

    def test_empty_header(self) -> None:
        with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=""):
            fields = parse_vcf_header("empty.vcf")
            assert fields == []

    def test_no_info_format_lines(self) -> None:
        header = "##fileformat=VCFv4.2\n#CHROM\tPOS\n"
        with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=header):
            fields = parse_vcf_header("minimal.vcf")
            assert fields == []

    def test_description_with_special_chars(self) -> None:
        header = "##INFO=<ID=ANN,Number=.,Type=String,Description=\"Func annotations: 'A | B'\">\n"
        with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=header):
            fields = parse_vcf_header("special.vcf")
            assert len(fields) == 1
            assert fields[0].description == "Func annotations: 'A | B'"


# ---------------------------------------------------------------------------
# format_fields_table
# ---------------------------------------------------------------------------
class TestFormatFieldsTable:
    def test_contains_info_and_format_sections(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields)
        assert "INFO fields (7 total):" in output
        assert "FORMAT fields (4 total):" in output

    def test_pattern_filter(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields, pattern="dbNSFP")
        assert "dbNSFP_REVEL_score" in output
        assert "GT" not in output

    def test_info_only(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields, info_only=True)
        assert "INFO fields" in output
        assert "FORMAT fields" not in output

    def test_format_only(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields, format_only=True)
        assert "FORMAT fields" in output
        assert "INFO fields" not in output

    def test_no_match(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields, pattern="nonexistent_xyz")
        assert output == "No matching fields found."

    def test_case_insensitive_filter(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_table(sample_fields, pattern="DBNSFP")
        assert "dbNSFP_REVEL_score" in output


# ---------------------------------------------------------------------------
# format_fields_json
# ---------------------------------------------------------------------------
class TestFormatFieldsJson:
    def test_valid_json(self, sample_fields: list[AnnotationField]) -> None:
        output = format_fields_json(sample_fields)
        data = json.loads(output)
        assert isinstance(data, list)
        assert len(data) == 11  # 7 INFO + 4 FORMAT

    def test_json_field_structure(self, sample_fields: list[AnnotationField]) -> None:
        data = json.loads(format_fields_json(sample_fields))
        first = data[0]
        assert "id" in first
        assert "source" in first
        assert "type" in first
        assert "number" in first
        assert "description" in first

    def test_json_with_filter(self, sample_fields: list[AnnotationField]) -> None:
        data = json.loads(format_fields_json(sample_fields, pattern="GT"))
        # Should match GT (FORMAT)
        assert any(d["id"] == "GT" for d in data)


# ---------------------------------------------------------------------------
# validate_requested_fields
# ---------------------------------------------------------------------------
class TestValidateRequestedFields:
    def test_all_valid_fixed_columns(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("CHROM POS REF ALT", sample_fields)
        assert missing == []

    def test_valid_info_fields(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("AC AF dbNSFP_REVEL_score", sample_fields)
        assert missing == []

    def test_ann_compound_field(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("ANN[0].GENE ANN[0].EFFECT", sample_fields)
        assert missing == []

    def test_gen_format_field(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("GEN[*].GT GEN[*].DP", sample_fields)
        assert missing == []

    def test_missing_info_field(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields(
            "CHROM dbNSFP_REVEL_Score",
            sample_fields,  # note: wrong case
        )
        assert "dbNSFP_REVEL_Score" in missing

    def test_missing_format_field(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("GEN[*].XYZ", sample_fields)
        assert "GEN[*].XYZ" in missing

    def test_gerp_special_chars(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("dbNSFP_GERP++_RS", sample_fields)
        assert missing == []

    def test_mixed_valid_and_invalid(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("CHROM AC nonexistent_field GEN[*].DP", sample_fields)
        assert missing == ["nonexistent_field"]

    def test_empty_fields(self, sample_fields: list[AnnotationField]) -> None:
        missing = validate_requested_fields("", sample_fields)
        assert missing == []

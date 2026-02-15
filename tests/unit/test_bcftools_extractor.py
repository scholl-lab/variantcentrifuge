"""
Unit tests for bcftools-based field extraction.

Tests the bcftools query extraction backend including:
- Format string building from config fields
- ANN subfield parsing from raw pipe-delimited strings
- NMD subfield parsing
- Missing value normalization
- Per-sample column naming
"""

import pandas as pd
import pytest

from variantcentrifuge.extractor import (
    build_bcftools_format_string,
    parse_ann_subfields,
    parse_nmd_subfields,
)


@pytest.mark.unit
class TestBcftoolsFormatString:
    """Test bcftools query format string generation."""

    def test_fixed_fields(self):
        """Test fixed VCF field mapping."""
        fields = ["CHROM", "POS", "REF", "ALT", "ID"]
        format_str, columns = build_bcftools_format_string(fields)

        assert "%CHROM" in format_str
        assert "%POS" in format_str
        assert "%REF" in format_str
        assert "%ALT" in format_str
        assert "%ID" in format_str

        assert columns == ["CHROM", "POS", "REF", "ALT", "ID"]

    def test_info_fields(self):
        """Test INFO field mapping."""
        fields = ["AC", "dbNSFP_REVEL_score", "ClinVar_CLNSIG"]
        format_str, columns = build_bcftools_format_string(fields)

        assert "%INFO/AC" in format_str
        assert "%INFO/dbNSFP_REVEL_score" in format_str
        assert "%INFO/ClinVar_CLNSIG" in format_str

        assert columns == ["AC", "dbNSFP_REVEL_score", "ClinVar_CLNSIG"]

    def test_ann_fields_single(self):
        """Test single ANN subfield extraction."""
        fields = ["CHROM", "POS", "ANN[0].GENE"]
        format_str, columns = build_bcftools_format_string(fields)

        # Should include %INFO/ANN once
        assert format_str.count("%INFO/ANN") == 1
        assert columns == ["CHROM", "POS", "ANN"]

    def test_ann_fields_multiple(self):
        """Test multiple ANN subfields - should deduplicate."""
        fields = [
            "CHROM",
            "POS",
            "ANN[0].GENE",
            "ANN[0].EFFECT",
            "ANN[0].IMPACT",
            "ANN[0].HGVS_C",
        ]
        format_str, columns = build_bcftools_format_string(fields)

        # Should include %INFO/ANN only once (deduplicated)
        assert format_str.count("%INFO/ANN") == 1
        assert columns == ["CHROM", "POS", "ANN"]

    def test_nmd_fields(self):
        """Test NMD field extraction."""
        fields = ["CHROM", "POS", "NMD[0].PERC"]
        format_str, columns = build_bcftools_format_string(fields)

        assert "%INFO/NMD" in format_str
        assert columns == ["CHROM", "POS", "NMD"]

    def test_per_sample_gt(self):
        """Test per-sample GT field extraction."""
        fields = ["CHROM", "POS", "GEN[*].GT"]
        samples = ["Sample1", "Sample2", "Sample3"]
        format_str, columns = build_bcftools_format_string(fields, samples)

        # Should have per-sample format string
        assert "[\\t%GT]" in format_str

        # Columns should include CHROM, POS, then per-sample GT
        assert columns[:2] == ["CHROM", "POS"]
        assert columns[2:] == ["GEN[0].GT", "GEN[1].GT", "GEN[2].GT"]

    def test_per_sample_multiple_fields(self):
        """Test multiple per-sample FORMAT fields."""
        fields = ["CHROM", "POS", "GEN[*].GT", "GEN[*].DP"]
        samples = ["Sample1", "Sample2"]
        format_str, columns = build_bcftools_format_string(fields, samples)

        # Should have both GT and DP in per-sample format
        assert "[\\t%GT\\t%DP]" in format_str

        # Columns should be interleaved: GT_0, DP_0, GT_1, DP_1
        # Actually they're in order: GEN[0].GT, GEN[0].DP, GEN[1].GT, GEN[1].DP
        assert "GEN[0].GT" in columns
        assert "GEN[0].DP" in columns
        assert "GEN[1].GT" in columns
        assert "GEN[1].DP" in columns

    def test_mixed_fields_config_json_example(self):
        """Test realistic field list from config.json."""
        fields = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "ID",
            "FILTER",
            "QUAL",
            "AC",
            "ANN[0].GENE",
            "ANN[0].FEATUREID",
            "ANN[0].EFFECT",
            "ANN[0].IMPACT",
            "ANN[0].HGVS_C",
            "ANN[0].HGVS_P",
            "NMD[0].PERC",
            "ANN[0].AA_POS",
            "ANN[0].AA_LEN",
            "dbNSFP_REVEL_score",
            "dbNSFP_CADD_phred",
            "splice_dbscSNV_rf_score",
            "ClinVar_CLNSIG",
            "GEN[*].GT",
        ]
        samples = ["Sample1", "Sample2"]
        format_str, columns = build_bcftools_format_string(fields, samples)

        # Fixed fields present
        assert "%CHROM" in format_str
        assert "%POS" in format_str

        # INFO fields present
        assert "%INFO/AC" in format_str
        assert "%INFO/dbNSFP_REVEL_score" in format_str

        # ANN only once
        assert format_str.count("%INFO/ANN") == 1

        # NMD present
        assert "%INFO/NMD" in format_str

        # Per-sample GT
        assert "[\\t%GT]" in format_str

        # Columns include raw ANN and NMD (to be parsed later)
        assert "ANN" in columns
        assert "NMD" in columns
        assert "GEN[0].GT" in columns
        assert "GEN[1].GT" in columns


@pytest.mark.unit
class TestAnnParsing:
    """Test ANN subfield parsing from raw pipe-delimited strings."""

    def test_parse_single_annotation(self):
        """Test parsing a single annotation with all fields."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["12345"],
                "ANN": [
                    "G|missense_variant|MODERATE|BRCA2|ENSG00000139618|transcript|"
                    "NM_000059.3|protein_coding|10/27|c.1234G>A|p.Arg412His|1234/5678|"
                    "1000/3000|412/1000|500|"
                ],
            }
        )

        fields = [
            "CHROM",
            "POS",
            "ANN[0].GENE",
            "ANN[0].EFFECT",
            "ANN[0].IMPACT",
            "ANN[0].HGVS_C",
            "ANN[0].HGVS_P",
            "ANN[0].AA_POS",
            "ANN[0].AA_LEN",
        ]

        result = parse_ann_subfields(df, fields)

        assert "ANN" not in result.columns  # Raw ANN should be removed
        assert result["ANN[0].GENE"].iloc[0] == "BRCA2"
        assert result["ANN[0].EFFECT"].iloc[0] == "missense_variant"
        assert result["ANN[0].IMPACT"].iloc[0] == "MODERATE"
        assert result["ANN[0].HGVS_C"].iloc[0] == "c.1234G>A"
        assert result["ANN[0].HGVS_P"].iloc[0] == "p.Arg412His"
        assert result["ANN[0].AA_POS"].iloc[0] == "412"
        assert result["ANN[0].AA_LEN"].iloc[0] == "1000"

    def test_parse_multiple_annotations_takes_first(self):
        """Test that multiple comma-separated annotations use only the first."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["12345"],
                "ANN": [
                    "G|missense_variant|MODERATE|BRCA2|||||||c.1234G>A|p.Arg412His|||500|,"
                    "G|synonymous_variant|LOW|TP53|||||||c.567C>T|p.Leu189=|||100|"
                ],
            }
        )

        fields = ["CHROM", "POS", "ANN[0].GENE", "ANN[0].EFFECT"]

        result = parse_ann_subfields(df, fields)

        # Should take first annotation only
        assert result["ANN[0].GENE"].iloc[0] == "BRCA2"
        assert result["ANN[0].EFFECT"].iloc[0] == "missense_variant"

    def test_parse_missing_ann(self):
        """Test handling of missing ANN field."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "ANN": [None]})

        fields = ["CHROM", "POS", "ANN[0].GENE", "ANN[0].EFFECT"]

        result = parse_ann_subfields(df, fields)

        # Missing values should be NA
        assert result["ANN[0].GENE"].iloc[0] == "NA"
        assert result["ANN[0].EFFECT"].iloc[0] == "NA"

    def test_parse_empty_ann(self):
        """Test handling of empty ANN string."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "ANN": [""]})

        fields = ["CHROM", "POS", "ANN[0].GENE", "ANN[0].EFFECT"]

        result = parse_ann_subfields(df, fields)

        # Empty should become NA
        assert result["ANN[0].GENE"].iloc[0] == "NA"
        assert result["ANN[0].EFFECT"].iloc[0] == "NA"

    def test_parse_malformed_ann_fewer_fields(self):
        """Test handling of ANN with fewer pipe fields than expected."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["12345"],
                "ANN": ["G|missense_variant|MODERATE"],  # Only 3 fields, not 16
            }
        )

        fields = ["CHROM", "POS", "ANN[0].GENE", "ANN[0].HGVS_C"]

        result = parse_ann_subfields(df, fields)

        # When split expands beyond available columns, pandas fills with NaN
        # Our fillna("NA") converts these to "NA"
        # Position 3 (GENE) is beyond the 3 fields present, so becomes NA
        assert result["ANN[0].GENE"].iloc[0] == "NA"  # Position 3 doesn't exist
        assert result["ANN[0].HGVS_C"].iloc[0] == "NA"  # Position 9 doesn't exist

    def test_parse_partial_fields(self):
        """Test parsing only some ANN subfields."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["12345"],
                "ANN": ["G|missense_variant|MODERATE|BRCA2|||||||c.1234G>A|p.Arg412His|||500|"],
            }
        )

        # Only request GENE and EFFECT
        fields = ["CHROM", "POS", "ANN[0].GENE", "ANN[0].EFFECT"]

        result = parse_ann_subfields(df, fields)

        # Should only have requested subfields (plus original columns)
        assert "ANN[0].GENE" in result.columns
        assert "ANN[0].EFFECT" in result.columns
        assert "ANN[0].HGVS_C" not in result.columns  # Not requested

    def test_no_ann_fields_requested(self):
        """Test when ANN is in DataFrame but no ANN subfields requested."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "ANN": ["some_annotation"]})

        fields = ["CHROM", "POS"]  # No ANN subfields

        result = parse_ann_subfields(df, fields)

        # ANN column should be removed
        assert "ANN" not in result.columns
        assert list(result.columns) == ["CHROM", "POS"]

    def test_no_ann_column(self):
        """Test when ANN column is not present."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"]})

        fields = ["CHROM", "POS", "ANN[0].GENE"]

        result = parse_ann_subfields(df, fields)

        # Should return DataFrame unchanged (no ANN column to parse)
        assert list(result.columns) == ["CHROM", "POS"]


@pytest.mark.unit
class TestNmdParsing:
    """Test NMD subfield parsing."""

    def test_parse_nmd_perc(self):
        """Test parsing NMD[0].PERC field."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "NMD": ["0.95|other|fields"]})

        fields = ["CHROM", "POS", "NMD[0].PERC"]

        result = parse_nmd_subfields(df, fields)

        assert "NMD" not in result.columns  # Raw NMD should be removed
        assert result["NMD[0].PERC"].iloc[0] == "0.95"

    def test_parse_missing_nmd(self):
        """Test handling of missing NMD field."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "NMD": [None]})

        fields = ["CHROM", "POS", "NMD[0].PERC"]

        result = parse_nmd_subfields(df, fields)

        # When NMD is None, fillna("") creates empty string, split returns empty string
        # This gets normalized to "NA" in the main extract function
        assert result["NMD[0].PERC"].iloc[0] == ""  # Empty string before normalization

    def test_no_nmd_column(self):
        """Test when NMD column is not present."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"]})

        fields = ["CHROM", "POS", "NMD[0].PERC"]

        result = parse_nmd_subfields(df, fields)

        # Should return DataFrame unchanged
        assert list(result.columns) == ["CHROM", "POS"]


@pytest.mark.unit
class TestMissingValueNormalization:
    """Test missing value normalization patterns."""

    def test_dot_to_na(self):
        """Test that '.' values are normalized to 'NA'."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["."], "ID": ["."]})

        # In the actual extraction, this happens via fillna and replace
        df = df.fillna("NA")
        df = df.replace(".", "NA")

        assert df["POS"].iloc[0] == "NA"
        assert df["ID"].iloc[0] == "NA"

    def test_empty_to_na(self):
        """Test that empty strings are normalized to 'NA'."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": ["12345"], "ID": [""]})

        df = df.replace("", "NA")

        assert df["ID"].iloc[0] == "NA"


@pytest.mark.unit
class TestColumnNaming:
    """Test per-sample column naming."""

    def test_per_sample_column_names(self):
        """Test that per-sample columns are named correctly."""
        fields = ["CHROM", "POS", "GEN[*].GT"]
        samples = ["Sample1", "Sample2", "Sample3"]
        _, columns = build_bcftools_format_string(fields, samples)

        # Should have CHROM, POS, then GEN[0].GT, GEN[1].GT, GEN[2].GT
        assert columns == ["CHROM", "POS", "GEN[0].GT", "GEN[1].GT", "GEN[2].GT"]

    def test_per_sample_no_samples_list(self):
        """Test per-sample fields without sample list."""
        fields = ["CHROM", "POS", "GEN[*].GT"]
        _, columns = build_bcftools_format_string(fields, None)

        # Should fall back to generic naming
        assert columns == ["CHROM", "POS", "GEN[*].GT"]

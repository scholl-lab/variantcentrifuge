"""
Comprehensive test suite for the unified annotation module.

Tests all aspects of the annotation system including:
- BED file annotation with interval trees
- Gene list annotation
- JSON gene data annotation
- Configuration validation
- DataFrame processing
- Error handling
"""

import pytest
import pandas as pd
import tempfile
import os
import json
from unittest.mock import patch
from variantcentrifuge.annotator import (
    load_custom_features,
    annotate_dataframe_with_features,
    validate_annotation_config,
    get_annotation_summary,
    _sanitize_name,
    _parse_bed_line,
    _load_bed_files,
    _load_gene_lists,
    _load_json_gene_data,
    _extract_genes_from_row,
    _find_gene_list_matches,
    _find_json_gene_matches,
    INTERVALTREE_AVAILABLE,
)


class TestSanitizeName:
    """Test name sanitization functionality."""

    def test_basic_sanitization(self):
        """Test basic name sanitization."""
        assert _sanitize_name("gene_list.txt") == "gene_list"
        assert _sanitize_name("/path/to/cancer-genes.txt") == "cancer_genes"
        assert _sanitize_name("my gene list.txt") == "my_gene_list"

    def test_special_characters(self):
        """Test handling of special characters."""
        assert _sanitize_name("genes@2024.csv") == "genes_2024"
        assert (
            _sanitize_name("list-with-dashes_and_underscores.txt")
            == "list_with_dashes_and_underscores"
        )

    def test_empty_or_invalid_names(self):
        """Test handling of empty or invalid names."""
        assert _sanitize_name("") == "unnamed"
        assert _sanitize_name("...") == "unnamed"
        assert _sanitize_name("___") == "unnamed"


class TestBedLineParsing:
    """Test BED file line parsing."""

    def test_valid_bed_lines(self):
        """Test parsing of valid BED lines."""
        # Standard 4-column BED
        result = _parse_bed_line("chr1\t1000\t2000\tpromoter_region")
        assert result == ("1", 1000, 2000, "promoter_region")

        # 3-column BED (minimum)
        result = _parse_bed_line("chr2\t5000\t6000")
        assert result == ("2", 5000, 6000, "region_5000_6000")

        # Remove chr prefix
        result = _parse_bed_line("chrX\t100\t200\tX_region")
        assert result == ("X", 100, 200, "X_region")

    def test_invalid_bed_lines(self):
        """Test handling of invalid BED lines."""
        assert _parse_bed_line("") is None
        assert _parse_bed_line("# comment") is None
        assert _parse_bed_line("track name=test") is None
        assert _parse_bed_line("chr1\t1000") is None  # Too few columns
        assert _parse_bed_line("chr1\tinvalid\t2000") is None  # Invalid coordinates


class TestBedFileLoading:
    """Test BED file loading functionality."""

    @pytest.fixture
    def temp_bed_file(self):
        """Create a temporary BED file for testing."""
        bed_content = """# Test BED file
chr1\t1000\t2000\tpromoter_A
chr1\t3000\t4000\tenhancer_B
chr2\t5000\t6000\tpromoter_C
chrX\t7000\t8000\tregion_X
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False) as f:
            f.write(bed_content)
            return f.name

    def test_load_bed_files_success(self, temp_bed_file):
        """Test successful BED file loading."""
        if not INTERVALTREE_AVAILABLE:
            pytest.skip("intervaltree not available")

        regions = _load_bed_files([temp_bed_file])

        # Check chromosomes are present
        assert "1" in regions
        assert "2" in regions
        assert "X" in regions

        # Check interval counts
        assert len(regions["1"]) == 2
        assert len(regions["2"]) == 1
        assert len(regions["X"]) == 1

        # Check specific intervals
        overlaps_1_1500 = list(regions["1"][1500])
        assert len(overlaps_1_1500) == 1
        assert "promoter_A" in overlaps_1_1500[0].data

        # Cleanup
        os.unlink(temp_bed_file)

    def test_load_bed_files_missing_file(self):
        """Test handling of missing BED files."""
        regions = _load_bed_files(["/nonexistent/file.bed"])
        assert regions == {}

    @patch("variantcentrifuge.annotator.INTERVALTREE_AVAILABLE", False)
    def test_load_bed_files_no_intervaltree(self):
        """Test handling when intervaltree is not available."""
        regions = _load_bed_files(["any_file.bed"])
        assert regions == {}


class TestGeneListLoading:
    """Test gene list loading functionality."""

    @pytest.fixture
    def temp_gene_list(self):
        """Create a temporary gene list file."""
        gene_content = """# Cancer genes
BRCA1
BRCA2
TP53
# Comment line
EGFR,KRAS,PIK3CA
"""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write(gene_content)
            return f.name

    def test_load_gene_lists_success(self, temp_gene_list):
        """Test successful gene list loading."""
        gene_lists = _load_gene_lists([temp_gene_list])

        list_name = _sanitize_name(temp_gene_list)
        assert list_name in gene_lists

        genes = gene_lists[list_name]
        expected_genes = {"BRCA1", "BRCA2", "TP53", "EGFR", "KRAS", "PIK3CA"}
        assert genes == expected_genes

        # Cleanup
        os.unlink(temp_gene_list)

    def test_load_gene_lists_missing_file(self):
        """Test handling of missing gene list files."""
        gene_lists = _load_gene_lists(["/nonexistent/genes.txt"])
        assert gene_lists == {}


class TestJsonGeneDataLoading:
    """Test JSON gene data loading functionality."""

    @pytest.fixture
    def temp_json_file(self):
        """Create a temporary JSON gene data file."""
        gene_data = [
            {
                "gene_symbol": "BRCA1",
                "panel": "Hereditary_Cancer",
                "inheritance": "Autosomal_Dominant",
                "pathogenicity": "High",
            },
            {
                "gene_symbol": "CFTR",
                "panel": "Cystic_Fibrosis",
                "inheritance": "Autosomal_Recessive",
                "pathogenicity": "High",
            },
        ]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(gene_data, f)
            return f.name

    def test_load_json_gene_data_success(self, temp_json_file):
        """Test successful JSON gene data loading."""
        mapping_config = (
            '{"identifier":"gene_symbol","dataFields":["panel","inheritance","pathogenicity"]}'
        )

        json_data = _load_json_gene_data([temp_json_file], mapping_config)

        assert "BRCA1" in json_data
        assert "CFTR" in json_data

        brca1_data = json_data["BRCA1"]
        assert brca1_data["panel"] == "Hereditary_Cancer"
        assert brca1_data["inheritance"] == "Autosomal_Dominant"
        assert brca1_data["pathogenicity"] == "High"

        # Cleanup
        os.unlink(temp_json_file)

    def test_load_json_gene_data_invalid_mapping(self, temp_json_file):
        """Test handling of invalid mapping configuration."""
        json_data = _load_json_gene_data([temp_json_file], "invalid json")
        assert json_data == {}

        # Cleanup
        os.unlink(temp_json_file)

    def test_load_json_gene_data_no_mapping(self, temp_json_file):
        """Test handling when no mapping is provided."""
        json_data = _load_json_gene_data([temp_json_file], "")
        assert json_data == {}

        # Cleanup
        os.unlink(temp_json_file)


class TestFeatureLoading:
    """Test the main feature loading function."""

    @pytest.fixture
    def comprehensive_config(self):
        """Create comprehensive test configuration with all annotation types."""
        # Create temporary files
        bed_content = "chr1\t1000\t2000\tpromoter\n"
        bed_file = tempfile.NamedTemporaryFile(mode="w", suffix=".bed", delete=False)
        bed_file.write(bed_content)
        bed_file.close()

        gene_content = "BRCA1\nTP53\n"
        gene_file = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
        gene_file.write(gene_content)
        gene_file.close()

        json_content = [{"gene_symbol": "EGFR", "panel": "Cancer"}]
        json_file = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(json_content, json_file)
        json_file.close()

        config = {
            "annotate_bed_files": [bed_file.name],
            "annotate_gene_lists": [gene_file.name],
            "annotate_json_genes": [json_file.name],
            "json_gene_mapping": '{"identifier":"gene_symbol","dataFields":["panel"]}',
        }

        return config, [bed_file.name, gene_file.name, json_file.name]

    def test_load_custom_features_comprehensive(self, comprehensive_config):
        """Test loading all types of custom features."""
        config, temp_files = comprehensive_config

        try:
            features = load_custom_features(config)

            # Check all feature types are loaded
            if INTERVALTREE_AVAILABLE:
                assert "1" in features["regions_by_chrom"]
            assert len(features["gene_lists"]) == 1
            assert "EGFR" in features["json_gene_data"]

        finally:
            # Cleanup
            for temp_file in temp_files:
                os.unlink(temp_file)

    def test_load_custom_features_empty_config(self):
        """Test loading with empty configuration."""
        features = load_custom_features({})

        assert features["regions_by_chrom"] == {}
        assert features["gene_lists"] == {}
        assert features["json_gene_data"] == {}


class TestGeneExtraction:
    """Test gene extraction from variant rows."""

    def test_extract_genes_single_gene(self):
        """Test extracting a single gene."""
        row = pd.Series({"GENE": "BRCA1", "other_col": "value"})
        genes = _extract_genes_from_row(row)
        assert genes == {"BRCA1"}

    def test_extract_genes_multiple_genes(self):
        """Test extracting multiple genes with different separators."""
        # Comma-separated
        row = pd.Series({"GENE": "BRCA1,BRCA2,TP53"})
        genes = _extract_genes_from_row(row)
        assert genes == {"BRCA1", "BRCA2", "TP53"}

        # Semicolon-separated
        row = pd.Series({"GENE": "EGFR;KRAS;PIK3CA"})
        genes = _extract_genes_from_row(row)
        assert genes == {"EGFR", "KRAS", "PIK3CA"}

    def test_extract_genes_empty_or_missing(self):
        """Test handling of empty or missing gene information."""
        row = pd.Series({"GENE": "", "other_col": "value"})
        genes = _extract_genes_from_row(row)
        assert genes == set()

        row = pd.Series({"other_col": "value"})
        genes = _extract_genes_from_row(row)
        assert genes == set()

    def test_extract_genes_alternative_columns(self):
        """Test extracting genes from alternative column names."""
        row = pd.Series({"Gene": "BRCA1", "gene_symbol": "TP53"})
        genes = _extract_genes_from_row(row)
        assert "BRCA1" in genes or "TP53" in genes  # Should find at least one


class TestAnnotationFunctions:
    """Test individual annotation functions."""

    def test_find_gene_list_matches(self):
        """Test gene list matching."""
        genes = {"BRCA1", "TP53"}
        gene_lists = {
            "cancer_genes": {"BRCA1", "BRCA2", "TP53"},
            "metabolic_genes": {"LDLR", "PCSK9"},
        }

        matches = _find_gene_list_matches(genes, gene_lists)
        assert "InGeneList=cancer_genes" in matches
        assert "InGeneList=metabolic_genes" not in matches

    def test_find_json_gene_matches(self):
        """Test JSON gene data matching."""
        genes = {"BRCA1", "CFTR"}
        json_data = {
            "BRCA1": {"panel": "Cancer", "inheritance": "AD"},
            "LDLR": {"panel": "Lipid", "inheritance": "AD"},
        }

        matches = _find_json_gene_matches(genes, json_data)
        assert "panel=Cancer" in matches
        assert "inheritance=AD" in matches
        assert "panel=Lipid" not in matches


class TestDataFrameAnnotation:
    """Test full DataFrame annotation functionality."""

    @pytest.fixture
    def sample_dataframe(self):
        """Create a sample DataFrame for testing."""
        data = {
            "CHROM": ["chr1", "chr2", "chr1"],
            "POS": [1500, 5500, 3500],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "A"],
            "GENE": ["BRCA1", "TP53,EGFR", "UNKNOWN"],
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def sample_features(self):
        """Create sample features for testing."""
        features = {
            "regions_by_chrom": {},
            "gene_lists": {"cancer_genes": {"BRCA1", "TP53", "EGFR"}},
            "json_gene_data": {
                "BRCA1": {"panel": "Hereditary_Cancer"},
                "TP53": {"panel": "Cancer_Predisposition"},
            },
        }

        # Add interval tree if available
        if INTERVALTREE_AVAILABLE:
            from intervaltree import IntervalTree

            features["regions_by_chrom"]["1"] = IntervalTree()
            features["regions_by_chrom"]["1"][1000:2000] = "promoter_region"

        return features

    def test_annotate_dataframe_comprehensive(self, sample_dataframe, sample_features):
        """Test comprehensive DataFrame annotation."""
        result_df = annotate_dataframe_with_features(sample_dataframe, sample_features)

        # Check that Custom_Annotation column was added
        assert "Custom_Annotation" in result_df.columns

        # Check first variant (BRCA1 at chr1:1500)
        row1_annotation = result_df.iloc[0]["Custom_Annotation"]
        assert "InGeneList=cancer_genes" in row1_annotation
        assert "panel=Hereditary_Cancer" in row1_annotation
        if INTERVALTREE_AVAILABLE:
            assert "Region=promoter_region" in row1_annotation

        # Check second variant (TP53,EGFR at chr2:5500)
        row2_annotation = result_df.iloc[1]["Custom_Annotation"]
        assert "InGeneList=cancer_genes" in row2_annotation
        assert "panel=Cancer_Predisposition" in row2_annotation

        # Check third variant (UNKNOWN gene)
        row3_annotation = result_df.iloc[2]["Custom_Annotation"]
        assert row3_annotation == ""  # No matches expected

    def test_annotate_dataframe_empty(self):
        """Test annotation of empty DataFrame."""
        empty_df = pd.DataFrame()
        features = {"regions_by_chrom": {}, "gene_lists": {}, "json_gene_data": {}}

        result_df = annotate_dataframe_with_features(empty_df, features)

        assert "Custom_Annotation" in result_df.columns
        assert len(result_df) == 0

    def test_annotate_dataframe_no_features(self, sample_dataframe):
        """Test annotation with no features."""
        features = {"regions_by_chrom": {}, "gene_lists": {}, "json_gene_data": {}}

        result_df = annotate_dataframe_with_features(sample_dataframe, features)

        assert "Custom_Annotation" in result_df.columns
        assert all(result_df["Custom_Annotation"] == "")


class TestConfigValidation:
    """Test configuration validation functionality."""

    def test_validate_annotation_config_valid(self):
        """Test validation of valid configuration."""
        # Create temporary files
        bed_file = tempfile.NamedTemporaryFile(delete=False)
        gene_file = tempfile.NamedTemporaryFile(delete=False)
        json_file = tempfile.NamedTemporaryFile(delete=False)

        try:
            config = {
                "annotate_bed_files": [bed_file.name],
                "annotate_gene_lists": [gene_file.name],
                "annotate_json_genes": [json_file.name],
                "json_gene_mapping": '{"identifier":"gene_symbol","dataFields":["panel"]}',
            }

            errors = validate_annotation_config(config)
            assert errors == []

        finally:
            os.unlink(bed_file.name)
            os.unlink(gene_file.name)
            os.unlink(json_file.name)

    def test_validate_annotation_config_missing_files(self):
        """Test validation with missing files."""
        config = {
            "annotate_bed_files": ["/nonexistent/file.bed"],
            "annotate_gene_lists": ["/nonexistent/genes.txt"],
            "annotate_json_genes": ["/nonexistent/data.json"],
            "json_gene_mapping": '{"identifier":"gene_symbol","dataFields":["panel"]}',
        }

        errors = validate_annotation_config(config)
        assert len(errors) == 3  # Three missing files
        assert all("not found" in error for error in errors)

    def test_validate_annotation_config_json_mapping_errors(self):
        """Test validation of JSON mapping errors."""
        # Missing mapping for JSON files
        config = {"annotate_json_genes": ["file.json"], "json_gene_mapping": ""}
        errors = validate_annotation_config(config)
        assert any("no json_gene_mapping provided" in error for error in errors)

        # Invalid JSON syntax
        config = {"annotate_json_genes": ["file.json"], "json_gene_mapping": "invalid json"}
        errors = validate_annotation_config(config)
        assert any("Invalid JSON syntax" in error for error in errors)

        # Missing identifier field
        config = {
            "annotate_json_genes": ["file.json"],
            "json_gene_mapping": '{"dataFields":["panel"]}',
        }
        errors = validate_annotation_config(config)
        assert any("must include 'identifier' field" in error for error in errors)


class TestAnnotationSummary:
    """Test annotation summary functionality."""

    def test_get_annotation_summary(self):
        """Test generation of annotation summary."""
        df = pd.DataFrame(
            {
                "CHROM": ["chr1", "chr2", "chr3"],
                "Custom_Annotation": ["Region=promoter;InGeneList=cancer", "panel=Cancer", ""],
            }
        )

        summary = get_annotation_summary(df)

        assert summary["total_variants"] == 3
        assert summary["annotated_variants"] == 2
        assert summary["annotation_rate"] == 2 / 3
        assert summary["annotation_types"]["Region"] == 1
        assert summary["annotation_types"]["InGeneList"] == 1
        assert summary["annotation_types"]["panel"] == 1

    def test_get_annotation_summary_missing_column(self):
        """Test summary generation when column is missing."""
        df = pd.DataFrame({"CHROM": ["chr1"]})

        summary = get_annotation_summary(df)
        assert "error" in summary


class TestIntegrationScenarios:
    """Test realistic integration scenarios."""

    def test_cancer_gene_panel_annotation(self):
        """Test annotation scenario for cancer gene panel analysis."""
        # Create test data
        df = pd.DataFrame(
            {
                "CHROM": ["chr17", "chr13", "chr3"],
                "POS": [43124096, 32398489, 179234297],
                "REF": ["A", "G", "C"],
                "ALT": ["T", "C", "T"],
                "GENE": ["BRCA1", "BRCA2", "PIK3CA"],
            }
        )

        # Create features
        features = {
            "regions_by_chrom": {},
            "gene_lists": {
                "hereditary_cancer": {"BRCA1", "BRCA2", "TP53", "PALB2"},
                "actionable_genes": {"PIK3CA", "EGFR", "KRAS"},
            },
            "json_gene_data": {
                "BRCA1": {"panel": "Hereditary_Breast_Cancer", "tier": "1"},
                "BRCA2": {"panel": "Hereditary_Breast_Cancer", "tier": "1"},
                "PIK3CA": {"panel": "Targeted_Therapy", "tier": "2"},
            },
        }

        # Apply annotation
        result_df = annotate_dataframe_with_features(df, features)

        # Verify results
        brca1_annotation = result_df.iloc[0]["Custom_Annotation"]
        assert "InGeneList=hereditary_cancer" in brca1_annotation
        assert "panel=Hereditary_Breast_Cancer" in brca1_annotation
        assert "tier=1" in brca1_annotation

        pik3ca_annotation = result_df.iloc[2]["Custom_Annotation"]
        assert "InGeneList=actionable_genes" in pik3ca_annotation
        assert "panel=Targeted_Therapy" in pik3ca_annotation
        assert "tier=2" in pik3ca_annotation


class TestJsonColumnAnnotation:
    """Test JSON annotation as separate columns functionality."""

    @pytest.fixture
    def sample_dataframe_for_columns(self):
        """Create a sample DataFrame for column-based testing."""
        data = {
            "CHROM": ["chr1", "chr2", "chr3", "chr4"],
            "POS": [1000, 2000, 3000, 4000],
            "REF": ["A", "G", "C", "T"],
            "ALT": ["T", "C", "A", "G"],
            "GENE": ["BRCA1", "TP53", "UNKNOWN", "EGFR"],
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def features_with_json_columns(self):
        """Create features with json_genes_as_columns enabled."""
        features = {
            "regions_by_chrom": {},
            "gene_lists": {"cancer_genes": {"BRCA1", "TP53"}},
            "json_gene_data": {
                "BRCA1": {"ngs": "panel1", "actionability": "Tier1", "inheritance": "AD"},
                "TP53": {"ngs": "panel2", "actionability": "Tier1", "inheritance": "AD"},
                "EGFR": {"ngs": "panel3", "actionability": "Tier2"},  # Missing inheritance
            },
            "json_genes_as_columns": True,
        }
        return features

    def test_json_as_columns_functionality(
        self, sample_dataframe_for_columns, features_with_json_columns
    ):
        """Test that JSON data is added as separate columns."""
        result_df = annotate_dataframe_with_features(
            sample_dataframe_for_columns, features_with_json_columns
        )

        # Check that new columns were added
        assert "ngs" in result_df.columns
        assert "actionability" in result_df.columns
        assert "inheritance" in result_df.columns

        # Check column values for BRCA1
        assert result_df.iloc[0]["ngs"] == "panel1"
        assert result_df.iloc[0]["actionability"] == "Tier1"
        assert result_df.iloc[0]["inheritance"] == "AD"

        # Check column values for TP53
        assert result_df.iloc[1]["ngs"] == "panel2"
        assert result_df.iloc[1]["actionability"] == "Tier1"
        assert result_df.iloc[1]["inheritance"] == "AD"

        # Check that UNKNOWN gene has empty/NA values
        assert pd.isna(result_df.iloc[2]["ngs"]) or result_df.iloc[2]["ngs"] == ""
        assert (
            pd.isna(result_df.iloc[2]["actionability"]) or result_df.iloc[2]["actionability"] == ""
        )

        # Check EGFR has inheritance as NA (missing field)
        assert result_df.iloc[3]["ngs"] == "panel3"
        assert result_df.iloc[3]["actionability"] == "Tier2"
        assert pd.isna(result_df.iloc[3]["inheritance"]) or result_df.iloc[3]["inheritance"] == ""

        # Verify that JSON data is NOT in Custom_Annotation
        for annotation in result_df["Custom_Annotation"]:
            assert "ngs=" not in annotation
            assert "actionability=" not in annotation
            assert "inheritance=" not in annotation

    def test_json_as_columns_with_other_annotations(
        self, sample_dataframe_for_columns, features_with_json_columns
    ):
        """Test that other annotations still work with column mode."""
        result_df = annotate_dataframe_with_features(
            sample_dataframe_for_columns, features_with_json_columns
        )

        # Check that gene list annotations still appear in Custom_Annotation
        assert "InGeneList=cancer_genes" in result_df.iloc[0]["Custom_Annotation"]  # BRCA1
        assert "InGeneList=cancer_genes" in result_df.iloc[1]["Custom_Annotation"]  # TP53
        assert result_df.iloc[2]["Custom_Annotation"] == ""  # UNKNOWN
        assert result_df.iloc[3]["Custom_Annotation"] == ""  # EGFR not in cancer_genes

    def test_json_default_behavior_regression(self, sample_dataframe_for_columns):
        """Test that default behavior (json_genes_as_columns=False) still works."""
        features = {
            "regions_by_chrom": {},
            "gene_lists": {},
            "json_gene_data": {
                "BRCA1": {"panel": "Cancer", "tier": "1"},
                "TP53": {"panel": "Cancer", "tier": "1"},
            },
            "json_genes_as_columns": False,  # Explicit false
        }

        result_df = annotate_dataframe_with_features(sample_dataframe_for_columns, features)

        # Check that JSON columns were NOT added
        assert "panel" not in result_df.columns
        assert "tier" not in result_df.columns

        # Check that JSON data IS in Custom_Annotation
        assert "panel=Cancer" in result_df.iloc[0]["Custom_Annotation"]
        assert "tier=1" in result_df.iloc[0]["Custom_Annotation"]
        assert "panel=Cancer" in result_df.iloc[1]["Custom_Annotation"]
        assert "tier=1" in result_df.iloc[1]["Custom_Annotation"]

    def test_json_as_columns_empty_data(self):
        """Test column mode with no JSON data."""
        df = pd.DataFrame({"CHROM": ["chr1"], "POS": [1000], "GENE": ["BRCA1"]})
        features = {
            "regions_by_chrom": {},
            "gene_lists": {},
            "json_gene_data": {},  # Empty
            "json_genes_as_columns": True,
        }

        result_df = annotate_dataframe_with_features(df, features)

        # Should not add any new columns beyond Custom_Annotation
        assert set(result_df.columns) == set(df.columns) | {"Custom_Annotation"}

    def test_json_as_columns_all_columns_present(self):
        """Test that all unique fields from all genes create columns."""
        df = pd.DataFrame(
            {"CHROM": ["chr1", "chr2"], "POS": [1000, 2000], "GENE": ["GENE1", "GENE2"]}
        )

        features = {
            "regions_by_chrom": {},
            "gene_lists": {},
            "json_gene_data": {
                "GENE1": {"field1": "value1", "field2": "value2"},
                "GENE2": {"field2": "value2b", "field3": "value3"},  # field1 missing, field3 new
            },
            "json_genes_as_columns": True,
        }

        result_df = annotate_dataframe_with_features(df, features)

        # All unique fields should be present as columns
        assert "field1" in result_df.columns
        assert "field2" in result_df.columns
        assert "field3" in result_df.columns

        # Check values
        assert result_df.iloc[0]["field1"] == "value1"
        assert result_df.iloc[0]["field2"] == "value2"
        assert pd.isna(result_df.iloc[0]["field3"]) or result_df.iloc[0]["field3"] == ""

        assert pd.isna(result_df.iloc[1]["field1"]) or result_df.iloc[1]["field1"] == ""
        assert result_df.iloc[1]["field2"] == "value2b"
        assert result_df.iloc[1]["field3"] == "value3"


if __name__ == "__main__":
    pytest.main([__file__])

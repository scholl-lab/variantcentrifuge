"""Tests for the main inheritance analyzer."""

import pytest
import pandas as pd
import json
from variantcentrifuge.inheritance.analyzer import (
    analyze_inheritance,
    prepare_variant_info,
    create_inheritance_details,
    get_inheritance_summary,
    filter_by_inheritance_pattern,
    export_inheritance_report,
    process_inheritance_output,
)
import tempfile
import os


class TestInheritanceAnalyzer:
    """Test the main inheritance analyzer."""

    @pytest.fixture
    def trio_pedigree(self):
        """Provide standard trio pedigree."""
        return {
            "father": {
                "family_id": "FAM1",
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "mother": {
                "family_id": "FAM1",
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "child": {
                "family_id": "FAM1",
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
        }

    @pytest.fixture
    def de_novo_variants_df(self):
        """Provide DataFrame with de novo variants."""
        data = {
            "CHROM": ["chr1", "chr2"],
            "POS": [1000, 2000],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE2"],
            "IMPACT": ["HIGH", "MODERATE"],
            "AF": [0.0001, 0.0005],
            "father": ["0/0", "0/0"],
            "mother": ["0/0", "0/0"],
            "child": ["0/1", "0/1"],
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def recessive_variants_df(self):
        """Provide DataFrame with recessive inheritance pattern."""
        data = {
            "CHROM": ["chr3"],
            "POS": [3000],
            "REF": ["C"],
            "ALT": ["T"],
            "GENE": ["GENE3"],
            "IMPACT": ["HIGH"],
            "AF": [0.001],
            "father": ["0/1"],  # Carrier
            "mother": ["0/1"],  # Carrier
            "child": ["1/1"],  # Affected
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def compound_het_variants_df(self):
        """Provide DataFrame with compound heterozygous variants."""
        data = {
            "CHROM": ["chr4", "chr4"],
            "POS": [4000, 5000],
            "REF": ["A", "G"],
            "ALT": ["C", "T"],
            "GENE": ["GENE4", "GENE4"],
            "IMPACT": ["HIGH", "HIGH"],
            "AF": [0.001, 0.002],
            "father": ["0/1", "0/0"],  # Has variant 1
            "mother": ["0/0", "0/1"],  # Has variant 2
            "child": ["0/1", "0/1"],  # Has both
        }
        return pd.DataFrame(data)

    def test_analyze_inheritance_de_novo(self, de_novo_variants_df, trio_pedigree):
        """Test inheritance analysis for de novo variants."""
        sample_list = ["father", "mother", "child"]

        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        # Check that columns were added
        assert "Inheritance_Pattern" in result_df.columns
        assert "Inheritance_Details" in result_df.columns

        # Check de novo detection
        assert all(result_df["Inheritance_Pattern"] == "de_novo")

        # Check details
        for _, row in result_df.iterrows():
            details = json.loads(row["Inheritance_Details"])
            assert details["primary_pattern"] == "de_novo"
            assert "de_novo" in details["all_patterns"]
            assert details["confidence"] > 0.5

    def test_analyze_inheritance_recessive(self, recessive_variants_df, trio_pedigree):
        """Test inheritance analysis for recessive pattern."""
        sample_list = ["father", "mother", "child"]

        result_df = analyze_inheritance(recessive_variants_df, trio_pedigree, sample_list)

        assert result_df.iloc[0]["Inheritance_Pattern"] == "autosomal_recessive"

        details = json.loads(result_df.iloc[0]["Inheritance_Details"])
        assert details["primary_pattern"] == "autosomal_recessive"
        assert details["affected_count"] == 1  # Only child affected

    def test_analyze_inheritance_compound_het(self, compound_het_variants_df, trio_pedigree):
        """Test inheritance analysis for compound heterozygous variants."""
        sample_list = ["father", "mother", "child"]

        result_df = analyze_inheritance(compound_het_variants_df, trio_pedigree, sample_list)

        # Both variants should be marked as compound het
        assert all(result_df["Inheritance_Pattern"] == "compound_heterozygous")

        # Check details for linkage
        for _, row in result_df.iterrows():
            details = json.loads(row["Inheritance_Details"])
            assert "compound_heterozygous" in details["all_patterns"]

            # Find child info
            child_info = next(
                s for s in details["samples_with_pattern"] if s["sample_id"] == "child"
            )
            assert "compound_het_partner" in child_info
            assert child_info["compound_het_gene"] == "GENE4"

    def test_analyze_inheritance_empty_df(self, trio_pedigree):
        """Test with empty DataFrame."""
        empty_df = pd.DataFrame()
        sample_list = ["father", "mother", "child"]

        result_df = analyze_inheritance(empty_df, trio_pedigree, sample_list)

        assert "Inheritance_Pattern" in result_df.columns
        assert "Inheritance_Details" in result_df.columns
        assert len(result_df) == 0

    def test_prepare_variant_info(self):
        """Test variant info preparation."""
        # Rare, deleterious variant
        row = pd.Series(
            {"AF": 0.0001, "gnomAD_AF": 0.00005, "IMPACT": "HIGH", "CLIN_SIG": "Pathogenic"}
        )

        info = prepare_variant_info(row)

        assert info["allele_frequency"] == 0.00005
        assert info["is_rare"] is True
        assert info["is_very_rare"] is True
        assert info["is_deleterious"] is True
        assert info["is_lof"] is True
        assert info["is_pathogenic"] is True

        # Common, benign variant
        row = pd.Series({"AF": 0.1, "IMPACT": "LOW", "CLIN_SIG": "Benign"})

        info = prepare_variant_info(row)

        assert info["allele_frequency"] == 0.1
        assert info["is_rare"] is False
        assert info["is_deleterious"] is False
        assert info["is_pathogenic"] is False

    def test_create_inheritance_details(self, trio_pedigree):
        """Test inheritance details creation."""
        row = pd.Series({"father": "0/1", "mother": "0/1", "child": "1/1"})

        all_patterns = ["autosomal_recessive"]
        best_pattern = "autosomal_recessive"
        confidence = 0.9
        comp_het_info = None
        sample_list = ["father", "mother", "child"]

        details = create_inheritance_details(
            row, best_pattern, all_patterns, confidence, comp_het_info, trio_pedigree, sample_list
        )

        assert details["primary_pattern"] == "autosomal_recessive"
        assert details["confidence"] == 0.9
        assert len(details["samples_with_pattern"]) == 3
        assert details["affected_count"] == 1
        assert details["carrier_count"] == 2

    def test_get_inheritance_summary(self, de_novo_variants_df, trio_pedigree):
        """Test summary generation."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        summary = get_inheritance_summary(result_df)

        assert summary["total_variants"] == 2
        assert summary["pattern_counts"]["de_novo"] == 2
        assert summary["de_novo_variants"] == 2

    def test_filter_by_inheritance_pattern(self, trio_pedigree):
        """Test filtering by inheritance pattern."""
        # Create mixed DataFrame
        data = {
            "CHROM": ["chr1", "chr2", "chr3"],
            "POS": [1000, 2000, 3000],
            "REF": ["A", "G", "C"],
            "ALT": ["T", "C", "G"],
            "GENE": ["GENE1", "GENE2", "GENE3"],
            "father": ["0/0", "0/1", "0/0"],
            "mother": ["0/0", "0/1", "0/0"],
            "child": ["0/1", "1/1", "0/0"],
        }
        df = pd.DataFrame(data)
        sample_list = ["father", "mother", "child"]

        result_df = analyze_inheritance(df, trio_pedigree, sample_list)

        # Filter for de novo only
        de_novo_df = filter_by_inheritance_pattern(result_df, ["de_novo"])
        assert len(de_novo_df) == 1
        assert de_novo_df.iloc[0]["POS"] == 1000

        # Filter for recessive patterns
        recessive_df = filter_by_inheritance_pattern(
            result_df, ["autosomal_recessive", "compound_heterozygous"]
        )
        assert len(recessive_df) == 1
        assert recessive_df.iloc[0]["POS"] == 2000

    def test_export_inheritance_report(self, de_novo_variants_df, trio_pedigree):
        """Test report export."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            export_inheritance_report(result_df, f.name)

            # Read back and verify
            with open(f.name, "r") as rf:
                report_data = json.load(rf)

            assert len(report_data) == 2
            assert all(item["inheritance_pattern"] == "de_novo" for item in report_data)

            # Clean up
            os.unlink(f.name)

    def test_missing_gene_column(self, trio_pedigree):
        """Test handling of missing GENE column."""
        data = {
            "CHROM": ["chr1"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "father": ["0/0"],
            "mother": ["0/0"],
            "child": ["0/1"],
        }
        df = pd.DataFrame(data)
        sample_list = ["father", "mother", "child"]

        # Should not crash
        result_df = analyze_inheritance(df, trio_pedigree, sample_list)
        assert "Inheritance_Pattern" in result_df.columns
        assert result_df.iloc[0]["Inheritance_Pattern"] == "de_novo"

    def test_missing_samples(self, de_novo_variants_df, trio_pedigree):
        """Test with missing samples in DataFrame."""
        # Remove mother column
        df = de_novo_variants_df.drop(columns=["mother"])
        sample_list = ["father", "child"]  # Mother not available

        # Should handle gracefully
        result_df = analyze_inheritance(df, trio_pedigree, sample_list)
        assert len(result_df) == 2
        # Patterns might be different without mother data
        assert "Inheritance_Pattern" in result_df.columns

    def test_process_inheritance_output_simple(self, de_novo_variants_df, trio_pedigree):
        """Test process_inheritance_output with simple mode."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        # Process with simple mode
        processed_df = process_inheritance_output(result_df, "simple")

        # Should have pattern but not details
        assert "Inheritance_Pattern" in processed_df.columns
        assert "Inheritance_Details" not in processed_df.columns
        assert len(processed_df) == len(result_df)

    def test_process_inheritance_output_full(self, de_novo_variants_df, trio_pedigree):
        """Test process_inheritance_output with full mode."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        # Process with full mode
        processed_df = process_inheritance_output(result_df, "full")

        # Should have both columns unchanged
        assert "Inheritance_Pattern" in processed_df.columns
        assert "Inheritance_Details" in processed_df.columns
        assert processed_df.equals(result_df)

    def test_process_inheritance_output_columns(self, de_novo_variants_df, trio_pedigree):
        """Test process_inheritance_output with columns mode."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(de_novo_variants_df, trio_pedigree, sample_list)

        # Process with columns mode
        processed_df = process_inheritance_output(result_df, "columns")

        # Should have pattern and new columns but not details
        assert "Inheritance_Pattern" in processed_df.columns
        assert "Inheritance_Details" not in processed_df.columns
        assert "Inheritance_Confidence" in processed_df.columns
        assert "Inheritance_Description" in processed_df.columns
        assert "Inheritance_Samples" in processed_df.columns

        # Check values
        for idx, row in processed_df.iterrows():
            assert row["Inheritance_Pattern"] == "de_novo"
            assert float(row["Inheritance_Confidence"]) > 0
            assert "de novo" in row["Inheritance_Description"].lower()
            assert "child(0/1)" in row["Inheritance_Samples"]

    def test_process_inheritance_output_compound_het_columns(
        self, compound_het_variants_df, trio_pedigree
    ):
        """Test process_inheritance_output columns mode with compound het variants."""
        sample_list = ["father", "mother", "child"]
        result_df = analyze_inheritance(compound_het_variants_df, trio_pedigree, sample_list)

        # Process with columns mode
        processed_df = process_inheritance_output(result_df, "columns")

        # Check that compound het partner info is included
        for idx, row in processed_df.iterrows():
            if "partner:" in row["Inheritance_Samples"]:
                # This is a compound het variant
                assert "compound" in row["Inheritance_Description"].lower()

    def test_process_inheritance_output_no_inheritance_columns(self):
        """Test process_inheritance_output with DataFrame missing inheritance columns."""
        # Create DataFrame without inheritance columns
        df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": [1000],
                "REF": ["A"],
                "ALT": ["T"],
            }
        )

        # Should return unchanged for all modes
        for mode in ["simple", "columns", "full"]:
            processed_df = process_inheritance_output(df, mode)
            assert processed_df.equals(df)

    def test_process_inheritance_output_invalid_json(self):
        """Test process_inheritance_output with invalid JSON in details."""
        # Create DataFrame with invalid JSON
        df = pd.DataFrame(
            {
                "Inheritance_Pattern": ["de_novo"],
                "Inheritance_Details": ["invalid json"],
            }
        )

        # Process with columns mode - should handle gracefully
        processed_df = process_inheritance_output(df, "columns")

        assert "Inheritance_Confidence" in processed_df.columns
        assert processed_df.iloc[0]["Inheritance_Confidence"] == ""
        assert processed_df.iloc[0]["Inheritance_Description"] == ""
        assert processed_df.iloc[0]["Inheritance_Samples"] == ""

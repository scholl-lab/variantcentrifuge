"""Integration tests for the complete 3-pass inheritance analysis workflow."""

import json

import pandas as pd

from variantcentrifuge.inheritance.analyzer import (
    analyze_inheritance,
    create_inheritance_details,
    filter_by_inheritance_pattern,
    get_inheritance_summary,
)


class TestFullInheritanceWorkflow:
    """Test the complete inheritance analysis workflow."""

    def test_single_sample_analysis(self):
        """Test inheritance analysis for single sample without pedigree."""
        # Create test data
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "Sample1": "0/1",
                    "IMPACT": "HIGH",
                    "AF": "0.001",
                },
                {
                    "CHROM": "2",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE2",
                    "Sample1": "1/1",
                    "IMPACT": "MODERATE",
                    "AF": "0.01",
                },
                {
                    "CHROM": "3",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE3",
                    "Sample1": "0/0",
                    "IMPACT": "LOW",
                    "AF": "0.1",
                },
            ]
        )

        # Empty pedigree for single sample
        pedigree_data = {
            "Sample1": {
                "sample_id": "Sample1",
                "father_id": "0",
                "mother_id": "0",
                "sex": "0",
                "affected_status": "2",
            }
        }
        sample_list = ["Sample1"]

        # Run analysis
        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Check results
        assert "Inheritance_Pattern" in result_df.columns
        assert "Inheritance_Details" in result_df.columns

        # Row 1: Het variant should be unknown
        assert result_df.iloc[0]["Inheritance_Pattern"] == "unknown"

        # Row 2: Homozygous variant
        assert result_df.iloc[1]["Inheritance_Pattern"] == "homozygous"

        # Row 3: Reference
        assert result_df.iloc[2]["Inheritance_Pattern"] == "reference"

        # Check details are valid JSON
        for _, row in result_df.iterrows():
            details = json.loads(row["Inheritance_Details"])
            assert "primary_pattern" in details
            assert "confidence" in details
            assert "samples_with_pattern" in details

    def test_compound_het_detection(self):
        """Test compound heterozygous detection in single gene."""
        # Create test data with multiple variants in same gene
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "Sample1": "0/1",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "Sample1": "0/1",
                },
                {
                    "CHROM": "1",
                    "POS": "3000",
                    "REF": "T",
                    "ALT": "A",
                    "GENE": "GENE1",
                    "Sample1": "0/1",
                },
                {
                    "CHROM": "2",
                    "POS": "4000",
                    "REF": "G",
                    "ALT": "C",
                    "GENE": "GENE2",
                    "Sample1": "0/1",
                },
            ]
        )

        pedigree_data = {}  # No pedigree data
        sample_list = ["Sample1"]

        # Run analysis
        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # GENE1 variants should be compound het possible
        gene1_patterns = result_df[result_df["GENE"] == "GENE1"]["Inheritance_Pattern"].tolist()
        assert all(p == "compound_heterozygous_possible_no_pedigree" for p in gene1_patterns)

        # GENE2 single variant should be unknown
        gene2_pattern = result_df[result_df["GENE"] == "GENE2"]["Inheritance_Pattern"].iloc[0]
        assert gene2_pattern == "unknown"

    def test_family_analysis_de_novo(self):
        """Test de novo detection in family."""
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/0",
                    "mother": "0/0",
                },
            ]
        )

        pedigree_data = {
            "child": {
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
        }
        sample_list = ["child", "father", "mother"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert result_df.iloc[0]["Inheritance_Pattern"] == "de_novo"

    def test_family_analysis_recessive(self):
        """Test recessive pattern detection in family."""
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "1/1",
                    "father": "0/1",
                    "mother": "0/1",
                },
            ]
        )

        pedigree_data = {
            "child": {
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
        }
        sample_list = ["child", "father", "mother"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert result_df.iloc[0]["Inheritance_Pattern"] == "autosomal_recessive"


class TestInheritanceSummary:
    """Test inheritance summary generation."""

    def test_get_inheritance_summary(self):
        """Test generating summary statistics."""
        df = pd.DataFrame(
            [
                {
                    "Inheritance_Pattern": "de_novo",
                    "Inheritance_Details": '{"confidence": 0.9, "samples_with_pattern": []}',
                },
                {
                    "Inheritance_Pattern": "compound_heterozygous",
                    "Inheritance_Details": '{"confidence": 0.85, "samples_with_pattern": '
                    '[{"compound_het_gene": "GENE1"}]}',
                },
                {
                    "Inheritance_Pattern": "unknown",
                    "Inheritance_Details": '{"confidence": 0.1, "samples_with_pattern": []}',
                },
            ]
        )

        summary = get_inheritance_summary(df)

        assert summary["total_variants"] == 3
        assert summary["pattern_counts"]["de_novo"] == 1
        assert summary["pattern_counts"]["compound_heterozygous"] == 1
        assert summary["pattern_counts"]["unknown"] == 1
        assert summary["high_confidence_patterns"] == 2  # de_novo and compound_het
        assert summary["de_novo_variants"] == 1
        assert "GENE1" in summary["compound_het_genes"]


class TestInheritanceFiltering:
    """Test filtering by inheritance patterns."""

    def test_filter_by_pattern(self):
        """Test filtering variants by inheritance pattern."""
        df = pd.DataFrame(
            [
                {
                    "variant": "var1",
                    "Inheritance_Pattern": "de_novo",
                    "Inheritance_Details": '{"confidence": 0.9}',
                },
                {
                    "variant": "var2",
                    "Inheritance_Pattern": "autosomal_recessive",
                    "Inheritance_Details": '{"confidence": 0.8}',
                },
                {
                    "variant": "var3",
                    "Inheritance_Pattern": "unknown",
                    "Inheritance_Details": '{"confidence": 0.1}',
                },
            ]
        )

        # Filter for de novo
        filtered = filter_by_inheritance_pattern(df, ["de_novo"])
        assert len(filtered) == 1
        assert filtered.iloc[0]["variant"] == "var1"

        # Filter for multiple patterns
        filtered = filter_by_inheritance_pattern(df, ["de_novo", "autosomal_recessive"])
        assert len(filtered) == 2

        # Filter with confidence threshold
        filtered = filter_by_inheritance_pattern(
            df, ["de_novo", "autosomal_recessive"], min_confidence=0.85
        )
        assert len(filtered) == 1  # Only de novo has confidence >= 0.85


class TestInheritanceDetails:
    """Test inheritance details creation."""

    def test_create_details_single_sample(self):
        """Test creating details for single sample."""
        row = pd.Series({"Sample1": "0/1", "CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T"})

        pedigree_data = {"Sample1": {"sample_id": "Sample1", "affected_status": "2"}}

        details = create_inheritance_details(
            row, "unknown", ["unknown"], 0.5, None, pedigree_data, ["Sample1"]
        )

        assert details["primary_pattern"] == "unknown"
        assert details["confidence"] == 0.5
        assert len(details["samples_with_pattern"]) == 1
        assert details["samples_with_pattern"][0]["sample_id"] == "Sample1"
        assert details["samples_with_pattern"][0]["affected"] is True
        assert details["affected_count"] == 1
        assert details["carrier_count"] == 0

    def test_create_details_with_compound_het(self):
        """Test creating details with compound het info."""
        row = pd.Series({"Sample1": "0/1", "CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T"})

        comp_het_info = {
            "Sample1": {
                "is_compound_het": True,
                "partner_variant": "1:2000:C>G",
                "gene": "GENE1",
                "inheritance_type": "trans",
            }
        }

        pedigree_data = {"Sample1": {"sample_id": "Sample1", "affected_status": "2"}}

        details = create_inheritance_details(
            row,
            "compound_heterozygous",
            ["compound_heterozygous"],
            0.9,
            comp_het_info,
            pedigree_data,
            ["Sample1"],
        )

        assert details["samples_with_pattern"][0]["compound_het_partner"] == "1:2000:C>G"
        assert details["samples_with_pattern"][0]["compound_het_gene"] == "GENE1"
        assert details["samples_with_pattern"][0]["compound_het_configuration"] == "trans"


class TestEdgeCases:
    """Test edge cases in inheritance analysis."""

    def test_empty_dataframe(self):
        """Test handling of empty DataFrame."""
        df = pd.DataFrame()
        pedigree_data = {}
        sample_list = []

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert len(result_df) == 0
        assert "Inheritance_Pattern" in result_df.columns
        assert "Inheritance_Details" in result_df.columns

    def test_missing_gene_column(self):
        """Test handling when GENE column is missing."""
        df = pd.DataFrame(
            [
                {"CHROM": "1", "POS": "1000", "REF": "A", "ALT": "T", "Sample1": "0/1"},
                {"CHROM": "1", "POS": "2000", "REF": "C", "ALT": "G", "Sample1": "0/1"},
            ]
        )

        pedigree_data = {}
        sample_list = ["Sample1"]

        # Should not crash, just skip compound het analysis
        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert all(result_df["Inheritance_Pattern"] == "unknown")

    def test_all_missing_genotypes(self):
        """Test when all genotypes are missing."""
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "Sample1": "./.",
                    "Sample2": "./.",
                },
            ]
        )

        pedigree_data = {
            "Sample1": {"sample_id": "Sample1", "affected_status": "2"},
            "Sample2": {"sample_id": "Sample2", "affected_status": "1"},
        }
        sample_list = ["Sample1", "Sample2"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert result_df.iloc[0]["Inheritance_Pattern"] == "unknown"

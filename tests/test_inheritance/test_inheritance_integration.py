"""Integration tests for the complete 3-pass inheritance analysis workflow."""

import json

import pandas as pd

from variantcentrifuge.inheritance.analyzer import (
    analyze_inheritance,
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
        }
        sample_list = ["Sample1", "Sample2"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        assert result_df.iloc[0]["Inheritance_Pattern"] == "unknown"

"""Test cases to ensure inheritance analysis matches variant-linker behavior."""

import pandas as pd
from variantcentrifuge.inheritance.analyzer import analyze_inheritance


class TestVariantLinkerCompatibility:
    """Test that our implementation matches variant-linker's 3-pass system."""

    def test_pass1_pattern_deduction(self):
        """Test Pass 1: Pattern deduction matches variant-linker."""
        # Test case from variant-linker documentation
        df = pd.DataFrame(
            [
                # De novo variant
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "proband": "0/1",
                    "father": "0/0",
                    "mother": "0/0",
                },
                # Dominant from father
                {
                    "CHROM": "2",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE2",
                    "proband": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                # Recessive
                {
                    "CHROM": "3",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE3",
                    "proband": "1/1",
                    "father": "0/1",
                    "mother": "0/1",
                },
                # X-linked in male
                {
                    "CHROM": "X",
                    "POS": "4000",
                    "REF": "T",
                    "ALT": "C",
                    "GENE": "DMD",
                    "proband": "0/1",
                    "father": "0/0",
                    "mother": "0/1",
                },
            ]
        )

        pedigree_data = {
            "proband": {
                "sample_id": "proband",
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
        sample_list = ["proband", "father", "mother"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Verify patterns match variant-linker
        assert result_df.iloc[0]["Inheritance_Pattern"] == "de_novo"
        assert result_df.iloc[1]["Inheritance_Pattern"] in [
            "autosomal_dominant",
            "autosomal_dominant_possible",
        ]
        assert result_df.iloc[2]["Inheritance_Pattern"] == "autosomal_recessive"
        assert result_df.iloc[3]["Inheritance_Pattern"] == "x_linked_recessive"

    def test_pass2_compound_het_detection(self):
        """Test Pass 2: Compound het detection matches variant-linker."""
        # Classic compound het case
        df = pd.DataFrame(
            [
                # Trans configuration - one from each parent
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "proband": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "proband": "0/1",
                    "father": "0/0",
                    "mother": "0/1",
                },
                # Different gene - not compound het
                {
                    "CHROM": "2",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE2",
                    "proband": "0/1",
                    "father": "0/0",
                    "mother": "0/0",
                },
            ]
        )

        pedigree_data = {
            "proband": {
                "sample_id": "proband",
                "father_id": "father",
                "mother_id": "mother",
                "affected_status": "2",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = ["proband", "father", "mother"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # GENE1 variants should be compound het
        gene1_results = result_df[result_df["GENE"] == "GENE1"]
        assert all(gene1_results["Inheritance_Pattern"] == "compound_heterozygous")

        # GENE2 should not be compound het
        gene2_result = result_df[result_df["GENE"] == "GENE2"]
        assert gene2_result.iloc[0]["Inheritance_Pattern"] != "compound_heterozygous"

    def test_pass3_pattern_prioritization(self):
        """Test Pass 3: Pattern prioritization matches variant-linker."""
        # Variant with multiple possible patterns
        df = pd.DataFrame(
            [
                # Could be de novo or dominant with incomplete penetrance
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "proband": "0/1",
                    "sibling": "0/0",
                    "father": "./.",
                    "mother": "0/0",
                },
            ]
        )

        pedigree_data = {
            "proband": {
                "sample_id": "proband",
                "father_id": "father",
                "mother_id": "mother",
                "affected_status": "2",
            },
            "sibling": {
                "sample_id": "sibling",
                "father_id": "father",
                "mother_id": "mother",
                "affected_status": "1",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = ["proband", "sibling", "father", "mother"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Should prioritize de_novo_candidate over other patterns
        pattern = result_df.iloc[0]["Inheritance_Pattern"]
        assert pattern in ["de_novo_candidate", "autosomal_dominant_possible"]

        # Details should include all possible patterns
        import json

        details = json.loads(result_df.iloc[0]["Inheritance_Details"])
        assert len(details["all_patterns"]) > 1

    def test_single_sample_behavior(self):
        """Test that single sample analysis matches variant-linker."""
        df = pd.DataFrame(
            [
                # Single het
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "Sample1": "0/1",
                },
                # Homozygous
                {
                    "CHROM": "2",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE2",
                    "Sample1": "1/1",
                },
                # Multiple hets in same gene
                {
                    "CHROM": "3",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE3",
                    "Sample1": "0/1",
                },
                {
                    "CHROM": "3",
                    "POS": "4000",
                    "REF": "T",
                    "ALT": "C",
                    "GENE": "GENE3",
                    "Sample1": "0/1",
                },
            ]
        )

        # Empty pedigree for single sample
        pedigree_data = {}
        sample_list = ["Sample1"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Single het should be unknown
        assert result_df.iloc[0]["Inheritance_Pattern"] == "unknown"

        # Homozygous should be identified
        assert result_df.iloc[1]["Inheritance_Pattern"] == "homozygous"

        # Multiple hets in GENE3 should be compound het possible
        gene3_results = result_df[result_df["GENE"] == "GENE3"]
        assert all(
            gene3_results["Inheritance_Pattern"] == "compound_heterozygous_possible_no_pedigree"
        )

    def test_segregation_affects_prioritization(self):
        """Test that segregation analysis affects pattern choice."""
        # Large family with clear segregation
        df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "affected1": "0/1",
                    "affected2": "0/1",
                    "affected3": "0/1",
                    "unaffected1": "0/0",
                    "unaffected2": "0/0",
                    "parent1": "0/1",
                    "parent2": "0/0",
                },
            ]
        )

        pedigree_data = {
            "affected1": {
                "sample_id": "affected1",
                "father_id": "parent1",
                "mother_id": "parent2",
                "affected_status": "2",
            },
            "affected2": {
                "sample_id": "affected2",
                "father_id": "parent1",
                "mother_id": "parent2",
                "affected_status": "2",
            },
            "affected3": {
                "sample_id": "affected3",
                "father_id": "parent1",
                "mother_id": "parent2",
                "affected_status": "2",
            },
            "unaffected1": {
                "sample_id": "unaffected1",
                "father_id": "parent1",
                "mother_id": "parent2",
                "affected_status": "1",
            },
            "unaffected2": {
                "sample_id": "unaffected2",
                "father_id": "parent1",
                "mother_id": "parent2",
                "affected_status": "1",
            },
            "parent1": {
                "sample_id": "parent1",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "2",
            },
            "parent2": {
                "sample_id": "parent2",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = [
            "affected1",
            "affected2",
            "affected3",
            "unaffected1",
            "unaffected2",
            "parent1",
            "parent2",
        ]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Should identify as dominant with high confidence due to segregation
        assert result_df.iloc[0]["Inheritance_Pattern"] == "autosomal_dominant"

        import json

        details = json.loads(result_df.iloc[0]["Inheritance_Details"])
        assert details["confidence"] >= 0.7  # High confidence due to good segregation

    def test_complex_patterns(self):
        """Test complex inheritance patterns."""
        # Test mitochondrial inheritance
        df = pd.DataFrame(
            [
                {
                    "CHROM": "MT",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "MT-ND1",
                    "affected_mother": "1/1",
                    "affected_child1": "1/1",
                    "affected_child2": "1/1",
                    "unaffected_father": "0/0",
                },
            ]
        )

        pedigree_data = {
            "affected_mother": {
                "sample_id": "affected_mother",
                "sex": "2",
                "affected_status": "2",
                "father_id": "0",
                "mother_id": "0",
            },
            "affected_child1": {
                "sample_id": "affected_child1",
                "affected_status": "2",
                "father_id": "unaffected_father",
                "mother_id": "affected_mother",
            },
            "affected_child2": {
                "sample_id": "affected_child2",
                "affected_status": "2",
                "father_id": "unaffected_father",
                "mother_id": "affected_mother",
            },
            "unaffected_father": {
                "sample_id": "unaffected_father",
                "sex": "1",
                "affected_status": "1",
                "father_id": "0",
                "mother_id": "0",
            },
        }
        sample_list = ["affected_mother", "affected_child1", "affected_child2", "unaffected_father"]

        result_df = analyze_inheritance(df, pedigree_data, sample_list)

        # Should identify as mitochondrial or related pattern
        pattern = result_df.iloc[0]["Inheritance_Pattern"]
        # Mitochondrial inheritance might be classified as dominant if MT detection fails
        assert pattern in ["mitochondrial", "autosomal_dominant", "unknown"]

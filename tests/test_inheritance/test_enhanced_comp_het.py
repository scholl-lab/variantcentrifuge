"""Test cases for enhanced compound heterozygous analysis."""

import pandas as pd
from variantcentrifuge.inheritance.comp_het import (
    analyze_gene_for_compound_het,
    determine_compound_het_type,
    create_variant_key,
)


class TestCompoundHetDetection:
    """Test compound heterozygous detection."""

    def setup_method(self):
        """Set up test data."""
        self.pedigree_data = {
            "child": {
                "sample_id": "child",
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
        self.sample_list = ["child", "father", "mother"]

    def test_confirmed_compound_het_trans(self):
        """Test confirmed compound het with clear trans configuration."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/0",
                    "mother": "0/1",
                },
            ]
        )

        results = analyze_gene_for_compound_het(gene_df, self.pedigree_data, self.sample_list)

        # Should find compound het for both variants
        assert len(results) == 2

        # Check first variant
        var1_key = "1:1000:A>T"
        assert var1_key in results
        assert "child" in results[var1_key]
        assert results[var1_key]["child"]["comp_het_type"] == "compound_heterozygous"
        assert results[var1_key]["child"]["inheritance_type"] == "trans"

    def test_compound_het_possible_no_pedigree(self):
        """Test compound het detection with no pedigree data."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "sample1": "0/1",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "sample1": "0/1",
                },
            ]
        )

        # Empty pedigree data
        pedigree_data = {}
        sample_list = ["sample1"]

        results = analyze_gene_for_compound_het(gene_df, pedigree_data, sample_list)

        # Should still identify possible compound het
        assert len(results) == 2
        var1_key = "1:1000:A>T"
        assert var1_key in results
        assert (
            results[var1_key]["sample1"]["comp_het_type"]
            == "compound_heterozygous_possible_no_pedigree"
        )

    def test_compound_het_possible_missing_genotypes(self):
        """Test compound het with missing parent genotypes."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "./.",
                    "mother": "0/1",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "./.",
                },
            ]
        )

        results = analyze_gene_for_compound_het(gene_df, self.pedigree_data, self.sample_list)

        # Should identify as possible due to missing data
        assert len(results) == 2
        var1_key = "1:1000:A>T"
        assert (
            results[var1_key]["child"]["comp_het_type"]
            == "compound_heterozygous_possible_missing_parent_genotypes"
        )

    def test_compound_het_ambiguous(self):
        """Test ambiguous compound het where both parents have both variants."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/1",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/1",
                },
            ]
        )

        results = analyze_gene_for_compound_het(gene_df, self.pedigree_data, self.sample_list)

        # Should identify as possible due to ambiguous inheritance
        assert len(results) == 2
        var1_key = "1:1000:A>T"
        assert results[var1_key]["child"]["comp_het_type"] == "compound_heterozygous_possible"
        assert results[var1_key]["child"]["inheritance_type"] == "ambiguous"

    def test_not_compound_het_cis(self):
        """Test cis configuration (both variants from same parent)."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",
                },
            ]
        )

        results = analyze_gene_for_compound_het(gene_df, self.pedigree_data, self.sample_list)

        # Should not identify as compound het (both from father)
        assert len(results) == 0 or (
            len(results) == 2
            and results.get("1:1000:A>T", {}).get("child", {}).get("comp_het_type")
            == "not_compound_heterozygous"
        )


class TestCompoundHetTypeDetection:
    """Test specific compound het type determination."""

    def test_determine_type_trans_confirmed(self):
        """Test determining confirmed trans configuration."""
        var1 = pd.Series({"father": "0/1", "mother": "0/0"})
        var2 = pd.Series({"father": "0/0", "mother": "0/1"})
        gene_df = pd.DataFrame([var1, var2])

        pedigree_data = {"child": {"father_id": "father", "mother_id": "mother"}}
        sample_list = ["child", "father", "mother"]

        trans_type, comp_het_type = determine_compound_het_type(
            "child", 0, 1, gene_df, pedigree_data, sample_list
        )

        assert trans_type == "trans"
        assert comp_het_type == "compound_heterozygous"

    def test_determine_type_cis(self):
        """Test determining cis configuration."""
        var1 = pd.Series({"father": "0/1", "mother": "0/0"})
        var2 = pd.Series({"father": "0/1", "mother": "0/0"})
        gene_df = pd.DataFrame([var1, var2])

        pedigree_data = {"child": {"father_id": "father", "mother_id": "mother"}}
        sample_list = ["child", "father", "mother"]

        trans_type, comp_het_type = determine_compound_het_type(
            "child", 0, 1, gene_df, pedigree_data, sample_list
        )

        assert trans_type == "cis"
        assert comp_het_type == "not_compound_heterozygous"

    def test_determine_type_no_pedigree(self):
        """Test type determination without pedigree."""
        var1 = pd.Series({})
        var2 = pd.Series({})
        gene_df = pd.DataFrame([var1, var2])

        pedigree_data = {}
        sample_list = ["sample1"]

        trans_type, comp_het_type = determine_compound_het_type(
            "sample1", 0, 1, gene_df, pedigree_data, sample_list
        )

        assert trans_type == "unknown"
        assert comp_het_type == "compound_heterozygous_possible_no_pedigree"


class TestVariantKey:
    """Test variant key generation."""

    def test_create_variant_key(self):
        """Test creating unique variant keys."""
        variant = pd.Series({"CHROM": "1", "POS": "12345", "REF": "A", "ALT": "T"})

        key = create_variant_key(variant)
        assert key == "1:12345:A>T"

    def test_create_variant_key_missing_data(self):
        """Test variant key with missing data."""
        variant = pd.Series({})

        key = create_variant_key(variant)
        assert key == "chr?:0:N>N"

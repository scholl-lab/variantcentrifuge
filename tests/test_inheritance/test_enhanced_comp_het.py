"""Test cases for enhanced compound heterozygous analysis."""

import numpy as np
import pandas as pd

from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized_vectorized,
    create_variant_key,
    find_potential_partners_vectorized,
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

        results = analyze_gene_for_compound_het_vectorized(gene_df, self.pedigree_data, self.sample_list)

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

        results = analyze_gene_for_compound_het_vectorized(gene_df, pedigree_data, sample_list)

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

        results = analyze_gene_for_compound_het_vectorized(gene_df, self.pedigree_data, self.sample_list)

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

        results = analyze_gene_for_compound_het_vectorized(gene_df, self.pedigree_data, self.sample_list)

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

        results = analyze_gene_for_compound_het_vectorized(gene_df, self.pedigree_data, self.sample_list)

        # Should not identify as compound het (both from father)
        assert len(results) == 0 or (
            len(results) == 2
            and results.get("1:1000:A>T", {}).get("child", {}).get("comp_het_type")
            == "not_compound_heterozygous"
        )


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


class TestPartnerBasedCompHet:
    """Test the new partner-based compound het detection logic."""

    def test_no_parents_three_variants(self):
        """Test Case 1: No parents, affected individual with 3 het variants."""
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
                {
                    "CHROM": "1",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE1",
                    "sample1": "0/1",
                },
            ]
        )

        # No pedigree data - sample will be treated as affected singleton
        pedigree_data = {}
        sample_list = ["sample1"]

        results = analyze_gene_for_compound_het_vectorized_vectorized(gene_df, pedigree_data, sample_list)

        # Each variant should list the other two as partners
        assert len(results) == 3

        # Check variant 1
        var1_key = "1:1000:A>T"
        assert var1_key in results
        partners_var1 = results[var1_key]["sample1"]["partner_variants"]
        assert len(partners_var1) == 2
        assert "1:2000:C>G" in partners_var1
        assert "1:3000:G>A" in partners_var1

        # Check variant 2
        var2_key = "1:2000:C>G"
        assert var2_key in results
        partners_var2 = results[var2_key]["sample1"]["partner_variants"]
        assert len(partners_var2) == 2
        assert "1:1000:A>T" in partners_var2
        assert "1:3000:G>A" in partners_var2

        # Check variant 3
        var3_key = "1:3000:G>A"
        assert var3_key in results
        partners_var3 = results[var3_key]["sample1"]["partner_variants"]
        assert len(partners_var3) == 2
        assert "1:1000:A>T" in partners_var3
        assert "1:2000:C>G" in partners_var3

    def test_both_parents_three_variants(self):
        """Test Case 2: Both parents present, child with 3 variants in trans config."""
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
                    "mother": "0/0",  # varA from father
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/0",
                    "mother": "0/1",  # varB from mother
                },
                {
                    "CHROM": "1",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",  # varC from father
                },
            ]
        )

        pedigree_data = {
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
        sample_list = ["child", "father", "mother"]

        results = analyze_gene_for_compound_het_vectorized_vectorized(gene_df, pedigree_data, sample_list)

        # Check partnerships based on trans configuration
        var1_key = "1:1000:A>T"  # from father
        var2_key = "1:2000:C>G"  # from mother
        var3_key = "1:3000:G>A"  # from father

        # varA (from father) should partner with varB (from mother)
        assert var1_key in results
        partners_var1 = results[var1_key]["child"]["partner_variants"]
        assert partners_var1 == ["1:2000:C>G"]

        # varB (from mother) should partner with both varA and varC (both from father)
        assert var2_key in results
        partners_var2 = results[var2_key]["child"]["partner_variants"]
        assert len(partners_var2) == 2
        assert "1:1000:A>T" in partners_var2
        assert "1:3000:G>A" in partners_var2

        # varC (from father) should partner with varB (from mother)
        assert var3_key in results
        partners_var3 = results[var3_key]["child"]["partner_variants"]
        assert partners_var3 == ["1:2000:C>G"]

    def test_one_parent_three_variants(self):
        """Test Case 3: Only mother present, child with 3 variants."""
        gene_df = pd.DataFrame(
            [
                {
                    "CHROM": "1",
                    "POS": "1000",
                    "REF": "A",
                    "ALT": "T",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "mother": "0/1",  # varA from mother
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "mother": "0/0",  # varB NOT from mother (from absent father)
                },
                {
                    "CHROM": "1",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "mother": "0/0",  # varC NOT from mother (from absent father)
                },
            ]
        )

        pedigree_data = {
            "child": {
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = ["child", "mother"]  # father not in sample list

        results = analyze_gene_for_compound_het_vectorized_vectorized(gene_df, pedigree_data, sample_list)

        # varA (from mother) should partner with varB and varC (not from mother)
        var1_key = "1:1000:A>T"
        assert var1_key in results
        partners_var1 = results[var1_key]["child"]["partner_variants"]
        assert len(partners_var1) == 2
        assert "1:2000:C>G" in partners_var1
        assert "1:3000:G>A" in partners_var1

        # varB (not from mother) should partner with varA (from mother)
        var2_key = "1:2000:C>G"
        assert var2_key in results
        partners_var2 = results[var2_key]["child"]["partner_variants"]
        assert partners_var2 == ["1:1000:A>T"]

        # varC (not from mother) should partner with varA (from mother)
        var3_key = "1:3000:G>A"
        assert var3_key in results
        partners_var3 = results[var3_key]["child"]["partner_variants"]
        assert partners_var3 == ["1:1000:A>T"]

    def test_find_potential_partners_vectorized(self):
        """Test the core partner-finding function directly."""
        # Test with clear trans configuration
        het_indices = np.array([0, 1, 2])
        father_genotypes = np.array([1, 0, 1])  # Has var 0 and 2
        mother_genotypes = np.array([0, 1, 0])  # Has var 1

        partners = find_potential_partners_vectorized(
            het_indices, father_genotypes, mother_genotypes
        )

        # Var 0 (from father) should partner with var 1 (from mother)
        assert 0 in partners
        assert partners[0] == [1]

        # Var 1 (from mother) should partner with vars 0 and 2 (from father)
        assert 1 in partners
        assert set(partners[1]) == {0, 2}

        # Var 2 (from father) should partner with var 1 (from mother)
        assert 2 in partners
        assert partners[2] == [1]

    def test_empty_partner_list(self):
        """Test case where no valid partners exist."""
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

        pedigree_data = {
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
        sample_list = ["child", "father", "mother"]

        results = analyze_gene_for_compound_het_vectorized_vectorized(gene_df, pedigree_data, sample_list)

        # Both variants from father - no trans configuration possible
        assert len(results) == 0

    def test_ambiguous_origin_variants(self):
        """Test variants with ambiguous origin (present in both parents)."""
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
                    "mother": "0/1",  # Present in both
                },
                {
                    "CHROM": "1",
                    "POS": "2000",
                    "REF": "C",
                    "ALT": "G",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/1",
                    "mother": "0/0",  # Clear: from father
                },
                {
                    "CHROM": "1",
                    "POS": "3000",
                    "REF": "G",
                    "ALT": "A",
                    "GENE": "GENE1",
                    "child": "0/1",
                    "father": "0/0",
                    "mother": "0/1",  # Clear: from mother
                },
            ]
        )

        pedigree_data = {
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
        sample_list = ["child", "father", "mother"]

        results = analyze_gene_for_compound_het_vectorized_vectorized(gene_df, pedigree_data, sample_list)

        # Ambiguous variant should partner with clear origin variants
        var1_key = "1:1000:A>T"
        assert var1_key in results
        partners_var1 = results[var1_key]["child"]["partner_variants"]
        assert len(partners_var1) == 2
        assert "1:2000:C>G" in partners_var1
        assert "1:3000:G>A" in partners_var1

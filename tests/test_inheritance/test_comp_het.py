"""Tests for compound heterozygous analyzer."""

import pytest
import pandas as pd
from variantcentrifuge.inheritance.comp_het import (
    analyze_gene_for_compound_het,
    find_compound_het_pairs,
    is_potential_compound_het,
    is_trans_configuration,
    create_variant_key,
    get_compound_het_summary,
)


class TestCompoundHetAnalyzer:

    @pytest.fixture
    def trio_pedigree(self):
        """Standard trio pedigree."""
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
    def gene_variants_comp_het(self):
        """DataFrame with compound het variants in a gene."""
        data = {
            "CHROM": ["chr1", "chr1"],
            "POS": [1000, 2000],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE1"],
            "father": ["0/1", "0/0"],  # Father has variant 1
            "mother": ["0/0", "0/1"],  # Mother has variant 2
            "child": ["0/1", "0/1"],  # Child has both
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def gene_variants_not_comp_het(self):
        """DataFrame with variants that are not compound het."""
        data = {
            "CHROM": ["chr1", "chr1"],
            "POS": [1000, 2000],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GENE": ["GENE1", "GENE1"],
            "father": ["0/1", "0/1"],  # Father has both variants
            "mother": ["0/0", "0/0"],  # Mother has neither
            "child": ["0/1", "0/1"],  # Child has both (cis configuration)
        }
        return pd.DataFrame(data)

    @pytest.fixture
    def gene_variants_complex(self):
        """DataFrame with multiple variants for complex testing."""
        data = {
            "CHROM": ["chr1", "chr1", "chr1", "chr1"],
            "POS": [1000, 2000, 3000, 4000],
            "REF": ["A", "G", "C", "T"],
            "ALT": ["T", "C", "G", "A"],
            "GENE": ["GENE1", "GENE1", "GENE1", "GENE1"],
            "father": ["0/1", "0/0", "0/1", "0/0"],
            "mother": ["0/0", "0/1", "0/0", "0/1"],
            "child": ["0/1", "0/1", "0/1", "0/1"],
        }
        return pd.DataFrame(data)

    def test_analyze_gene_for_compound_het(self, gene_variants_comp_het, trio_pedigree):
        """Test basic compound het analysis."""
        sample_list = ["father", "mother", "child"]

        results = analyze_gene_for_compound_het(gene_variants_comp_het, trio_pedigree, sample_list)

        # Should identify both variants as part of compound het
        assert len(results) == 2

        # Check variant keys
        var1_key = "chr1:1000:A>T"
        var2_key = "chr1:2000:G>C"

        assert var1_key in results
        assert var2_key in results

        # Check child has compound het for both variants
        assert "child" in results[var1_key]
        assert results[var1_key]["child"]["is_compound_het"] is True
        assert results[var1_key]["child"]["partner_variant"] == var2_key
        assert results[var1_key]["child"]["inheritance_type"] == "trans"

        assert "child" in results[var2_key]
        assert results[var2_key]["child"]["is_compound_het"] is True
        assert results[var2_key]["child"]["partner_variant"] == var1_key

    def test_find_compound_het_pairs(self, gene_variants_comp_het, trio_pedigree):
        """Test finding compound het pairs."""
        sample_list = ["father", "mother", "child"]

        pairs = find_compound_het_pairs("child", gene_variants_comp_het, trio_pedigree, sample_list)

        assert len(pairs) == 1
        assert pairs[0] == (0, 1) or pairs[0] == (1, 0)

    def test_is_potential_compound_het(self, gene_variants_comp_het, trio_pedigree):
        """Test potential compound het detection."""
        sample_list = ["father", "mother", "child"]

        # Test true compound het
        is_comp_het = is_potential_compound_het(
            "child", 0, 1, gene_variants_comp_het, trio_pedigree, sample_list
        )
        assert is_comp_het is True

        # Test with parents (should not be compound het)
        is_comp_het = is_potential_compound_het(
            "father", 0, 1, gene_variants_comp_het, trio_pedigree, sample_list
        )
        assert is_comp_het is False

    def test_is_trans_configuration(self, gene_variants_comp_het, trio_pedigree):
        """Test trans configuration detection."""
        # True trans configuration
        is_trans = is_trans_configuration("child", 0, 1, gene_variants_comp_het, trio_pedigree)
        assert is_trans is True

        # Test with cis configuration
        gene_variants_cis = gene_variants_comp_het.copy()
        gene_variants_cis.loc[0, "father"] = "0/1"
        gene_variants_cis.loc[0, "mother"] = "0/1"
        gene_variants_cis.loc[1, "father"] = "0/1"
        gene_variants_cis.loc[1, "mother"] = "0/1"

        is_trans = is_trans_configuration("child", 0, 1, gene_variants_cis, trio_pedigree)
        assert is_trans is False

    def test_no_compound_het(self, gene_variants_not_comp_het, trio_pedigree):
        """Test when variants are not compound het."""
        sample_list = ["father", "mother", "child"]

        results = analyze_gene_for_compound_het(
            gene_variants_not_comp_het, trio_pedigree, sample_list
        )

        # Should not identify compound het for child
        if results:
            for var_key in results:
                if "child" in results[var_key]:
                    assert results[var_key]["child"]["inheritance_type"] == "unknown"

    def test_multiple_compound_het_pairs(self, gene_variants_complex, trio_pedigree):
        """Test with multiple potential compound het pairs."""
        sample_list = ["father", "mother", "child"]

        pairs = find_compound_het_pairs("child", gene_variants_complex, trio_pedigree, sample_list)

        # Should find multiple pairs
        assert len(pairs) > 1

        # Verify pairs are valid
        for var1_idx, var2_idx in pairs:
            assert var1_idx != var2_idx
            assert var1_idx >= 0 and var1_idx < len(gene_variants_complex)
            assert var2_idx >= 0 and var2_idx < len(gene_variants_complex)

    def test_create_variant_key(self):
        """Test variant key creation."""
        variant = pd.Series({"CHROM": "chr1", "POS": 12345, "REF": "A", "ALT": "T"})

        key = create_variant_key(variant)
        assert key == "chr1:12345:A>T"

        # Test with missing values
        variant_missing = pd.Series({})
        key = create_variant_key(variant_missing)
        assert key == "chr?:0:N>N"

    def test_get_compound_het_summary(self, gene_variants_comp_het, trio_pedigree):
        """Test compound het summary generation."""
        sample_list = ["father", "mother", "child"]

        results = analyze_gene_for_compound_het(gene_variants_comp_het, trio_pedigree, sample_list)

        summary = get_compound_het_summary(results, "child")

        assert summary["total_comp_het_variants"] == 2
        assert len(summary["comp_het_pairs"]) == 1
        assert "GENE1" in summary["genes_with_comp_het"]

        # Check pair details
        pair = summary["comp_het_pairs"][0]
        assert pair["gene"] == "GENE1"
        assert pair["configuration"] == "trans"

    def test_single_variant_gene(self, trio_pedigree):
        """Test gene with only one variant (no compound het possible)."""
        data = {
            "CHROM": ["chr1"],
            "POS": [1000],
            "REF": ["A"],
            "ALT": ["T"],
            "GENE": ["GENE1"],
            "father": ["0/1"],
            "mother": ["0/0"],
            "child": ["0/1"],
        }
        gene_df = pd.DataFrame(data)
        sample_list = ["father", "mother", "child"]

        results = analyze_gene_for_compound_het(gene_df, trio_pedigree, sample_list)

        # Should return empty results
        assert results == {}

    def test_missing_parent_data(self, gene_variants_comp_het):
        """Test compound het analysis with missing parent data."""
        # Pedigree with no parent info
        orphan_pedigree = {
            "child": {
                "family_id": "FAM1",
                "sample_id": "child",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "2",
            }
        }

        sample_list = ["child"]

        results = analyze_gene_for_compound_het(
            gene_variants_comp_het[["CHROM", "POS", "REF", "ALT", "GENE", "child"]],
            orphan_pedigree,
            sample_list,
        )

        # Should still identify potential compound hets if child is affected
        assert len(results) > 0

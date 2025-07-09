"""Tests for genotype utility functions."""

from variantcentrifuge.genotype_utils import (
    could_be_de_novo,
    get_allele_count,
    get_genotype_type,
    is_het,
    is_hom_alt,
    is_mendelian_consistent,
    is_missing,
    is_phased,
    is_ref,
    is_variant,
    merge_genotypes,
    parse_genotype,
)


class TestGenotypeUtils:
    """Test genotype utility functions."""

    def test_parse_genotype(self):
        """Test parsing various genotype formats."""
        assert parse_genotype("0/0") == (0, 0)
        assert parse_genotype("0/1") == (0, 1)
        assert parse_genotype("1/1") == (1, 1)
        assert parse_genotype("1|0") == (1, 0)
        assert parse_genotype("./.") == (None, None)
        assert parse_genotype("0/.") == (0, None)
        assert parse_genotype("./1") == (None, 1)
        assert parse_genotype("") == (None, None)
        assert parse_genotype("invalid") == (None, None)
        assert parse_genotype("2/3") == (2, 3)  # Multi-allelic

    def test_is_het(self):
        """Test heterozygous genotype detection."""
        assert is_het("0/1") is True
        assert is_het("1/0") is True
        assert is_het("0|1") is True
        assert is_het("1|0") is True
        assert is_het("0/0") is False
        assert is_het("1/1") is False
        assert is_het("./.") is False
        assert is_het("0/.") is False
        assert is_het("2/3") is False  # Multi-allelic het

    def test_is_hom_alt(self):
        """Test homozygous alternate detection."""
        assert is_hom_alt("1/1") is True
        assert is_hom_alt("1|1") is True
        assert is_hom_alt("0/0") is False
        assert is_hom_alt("0/1") is False
        assert is_hom_alt("./.") is False
        assert is_hom_alt("2/2") is False  # Multi-allelic

    def test_is_ref(self):
        """Test homozygous reference detection."""
        assert is_ref("0/0") is True
        assert is_ref("0|0") is True
        assert is_ref("0/1") is False
        assert is_ref("1/1") is False
        assert is_ref("./.") is False
        assert is_ref("0/.") is False

    def test_is_variant(self):
        """Test variant detection (any alt allele)."""
        assert is_variant("0/1") is True
        assert is_variant("1/0") is True
        assert is_variant("1/1") is True
        assert is_variant("0/0") is False
        assert is_variant("./.") is False
        assert is_variant("0/.") is False

    def test_is_missing(self):
        """Test missing genotype detection."""
        assert is_missing("./.") is True
        assert is_missing("0/.") is True
        assert is_missing("./1") is True
        assert is_missing("") is True
        assert is_missing(None) is True
        assert is_missing("0/0") is False
        assert is_missing("0/1") is False
        assert is_missing("1/1") is False

    def test_get_allele_count(self):
        """Test alternate allele counting."""
        assert get_allele_count("0/0") == 0
        assert get_allele_count("0/1") == 1
        assert get_allele_count("1/0") == 1
        assert get_allele_count("1/1") == 2
        assert get_allele_count("./.") == 0
        assert get_allele_count("0/.") == 0
        assert get_allele_count("./1") == 0

    def test_is_phased(self):
        """Test phased genotype detection."""
        assert is_phased("0|1") is True
        assert is_phased("1|0") is True
        assert is_phased("1|1") is True
        assert is_phased("0/1") is False
        assert is_phased("./.") is False
        assert is_phased("") is False

    def test_get_genotype_type(self):
        """Test genotype type classification."""
        assert get_genotype_type("0/0") == "ref"
        assert get_genotype_type("0/1") == "het"
        assert get_genotype_type("1/0") == "het"
        assert get_genotype_type("1/1") == "hom_alt"
        assert get_genotype_type("./.") == "missing"
        assert get_genotype_type("0/.") == "missing"
        assert get_genotype_type("") == "missing"

    def test_is_mendelian_consistent(self):
        """Test Mendelian consistency checking."""
        # Consistent cases
        assert is_mendelian_consistent("0/0", "0/0", "0/0") is True
        assert is_mendelian_consistent("0/1", "0/0", "0/1") is True
        assert is_mendelian_consistent("0/1", "0/1", "0/0") is True
        assert is_mendelian_consistent("0/1", "0/1", "0/1") is True
        assert is_mendelian_consistent("1/1", "0/1", "0/1") is True
        assert is_mendelian_consistent("1/1", "1/1", "1/1") is True

        # Inconsistent cases
        assert is_mendelian_consistent("1/1", "0/0", "0/0") is False
        assert is_mendelian_consistent("0/1", "0/0", "0/0") is False

        # Missing data - assume consistent
        assert is_mendelian_consistent("0/1", "./.", "0/0") is True
        assert is_mendelian_consistent("./.", "0/0", "0/0") is True

    def test_could_be_de_novo(self):
        """Test de novo variant detection."""
        # True de novo cases
        assert could_be_de_novo("0/1", "0/0", "0/0") is True
        assert could_be_de_novo("1/1", "0/0", "0/0") is True

        # Not de novo
        assert could_be_de_novo("0/0", "0/0", "0/0") is False
        assert could_be_de_novo("0/1", "0/1", "0/0") is False
        assert could_be_de_novo("0/1", "0/0", "0/1") is False
        assert could_be_de_novo("1/1", "0/1", "0/1") is False

        # Missing data
        assert could_be_de_novo("0/1", "./.", "0/0") is False
        assert could_be_de_novo("0/1", "0/0", "./.") is False

    def test_merge_genotypes(self):
        """Test merging multiple genotypes."""
        # All same
        assert merge_genotypes(["0/0", "0/0", "0/0"]) == "0/0"
        assert merge_genotypes(["1/1", "1/1", "1/1"]) == "1/1"

        # Mixed
        assert merge_genotypes(["0/0", "0/1", "0/1"]) == "0/1"
        assert merge_genotypes(["0/1", "1/1"]) == "1/1" or merge_genotypes(["0/1", "1/1"]) == "0/1"

        # With missing
        assert merge_genotypes(["./.", "./.", "./."]) == "./."
        assert merge_genotypes(["0/0", "./.", "0/0"]) == "0/0"
        assert merge_genotypes(["0/1", "./."]) == "0/1"

        # Empty
        assert merge_genotypes([]) == "./."

        # Single genotype
        assert merge_genotypes(["0/1"]) == "0/1"

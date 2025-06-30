"""Tests for inheritance pattern deducer."""

import pytest
from variantcentrifuge.inheritance.deducer import (
    deduce_patterns_for_variant,
    deduce_patterns_for_sample,
    is_dominant_pattern,
    is_recessive_pattern,
    is_x_linked_pattern,
    check_segregation,
    get_inheritance_info,
)


class TestPatternDeducer:

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
    def extended_pedigree(self):
        """Extended family pedigree."""
        return {
            "grandfather": {
                "family_id": "FAM1",
                "sample_id": "grandfather",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "2",
            },
            "grandmother": {
                "family_id": "FAM1",
                "sample_id": "grandmother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "father": {
                "family_id": "FAM1",
                "sample_id": "father",
                "father_id": "grandfather",
                "mother_id": "grandmother",
                "sex": "1",
                "affected_status": "2",
            },
            "mother": {
                "family_id": "FAM1",
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
            "child1": {
                "family_id": "FAM1",
                "sample_id": "child1",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "1",
                "affected_status": "2",
            },
            "child2": {
                "family_id": "FAM1",
                "sample_id": "child2",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "2",
                "affected_status": "1",
            },
        }

    def test_de_novo_pattern(self, trio_pedigree):
        """Test de novo variant detection."""
        # De novo variant - child has variant, parents don't
        variant_row = {"father": "0/0", "mother": "0/0", "child": "0/1"}
        sample_list = ["father", "mother", "child"]

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        assert "de_novo" in patterns

        # Specifically for the child
        child_patterns = deduce_patterns_for_sample(
            "child", variant_row, trio_pedigree, sample_list
        )
        assert "de_novo" in child_patterns

    def test_autosomal_dominant_pattern(self, extended_pedigree):
        """Test autosomal dominant pattern detection."""
        # Dominant pattern - affected individuals have variant
        variant_row = {
            "grandfather": "0/1",  # Affected, has variant
            "grandmother": "0/0",  # Unaffected, no variant
            "father": "0/1",  # Affected, has variant
            "mother": "0/0",  # Unaffected, no variant
            "child1": "0/1",  # Affected, has variant
            "child2": "0/0",  # Unaffected, no variant
        }
        sample_list = list(extended_pedigree.keys())

        # Check father
        assert is_dominant_pattern("father", variant_row, extended_pedigree, sample_list)

        # Check child1
        assert is_dominant_pattern("child1", variant_row, extended_pedigree, sample_list)

        patterns = deduce_patterns_for_variant(variant_row, extended_pedigree, sample_list)
        assert "autosomal_dominant" in patterns

    def test_autosomal_recessive_pattern(self, trio_pedigree):
        """Test autosomal recessive pattern detection."""
        # Recessive pattern - affected child is homozygous, parents are carriers
        variant_row = {
            "father": "0/1",  # Carrier
            "mother": "0/1",  # Carrier
            "child": "1/1",  # Affected, homozygous
        }
        sample_list = ["father", "mother", "child"]

        assert is_recessive_pattern("child", variant_row, trio_pedigree, sample_list)

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        assert "autosomal_recessive" in patterns

    def test_x_linked_pattern(self, trio_pedigree):
        """Test X-linked pattern detection."""
        # X-linked pattern - affected male child, carrier mother
        variant_row = {
            "father": "0/0",  # Unaffected father
            "mother": "0/1",  # Carrier mother
            "child": "0/1",  # Affected male child (hemizygous)
        }
        sample_list = ["father", "mother", "child"]

        assert is_x_linked_pattern("child", variant_row, trio_pedigree, sample_list)

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        assert "x_linked_recessive" in patterns

    def test_no_pattern(self, trio_pedigree):
        """Test when no specific pattern is found."""
        # All individuals are reference
        variant_row = {"father": "0/0", "mother": "0/0", "child": "0/0"}
        sample_list = ["father", "mother", "child"]

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        assert patterns == []

    def test_unknown_pattern(self, trio_pedigree):
        """Test unknown pattern when variant present but no clear inheritance."""
        # Unaffected parent has variant
        trio_pedigree["father"]["affected_status"] = "1"  # Unaffected
        variant_row = {
            "father": "0/1",  # Unaffected with variant
            "mother": "0/0",
            "child": "0/1",  # Affected with variant
        }
        sample_list = ["father", "mother", "child"]

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        assert "unknown" in patterns or len(patterns) > 0

    def test_check_segregation(self, extended_pedigree):
        """Test segregation checking."""
        # Perfect dominant segregation
        variant_row = {
            "grandfather": "0/1",  # Affected
            "grandmother": "0/0",  # Unaffected
            "father": "0/1",  # Affected
            "mother": "0/0",  # Unaffected
            "child1": "0/1",  # Affected
            "child2": "0/0",  # Unaffected
        }
        sample_list = list(extended_pedigree.keys())

        assert check_segregation(variant_row, extended_pedigree, sample_list, "autosomal_dominant")

    def test_get_inheritance_info(self, trio_pedigree):
        """Test getting detailed inheritance information."""
        variant_row = {"father": "0/1", "mother": "0/1", "child": "1/1"}
        patterns = ["autosomal_recessive"]

        info = get_inheritance_info("child", patterns, variant_row, trio_pedigree)

        assert info["sample_id"] == "child"
        assert info["patterns"] == ["autosomal_recessive"]
        assert info["genotype"] == "1/1"
        assert info["affected_status"] is True
        assert info["father_genotype"] == "0/1"
        assert info["mother_genotype"] == "0/1"
        assert info["father_affected"] is False
        assert info["mother_affected"] is False

    def test_missing_samples(self, trio_pedigree):
        """Test handling of missing samples in variant data."""
        # Mother not in variant data
        variant_row = {"father": "0/0", "child": "0/1"}
        sample_list = ["father", "child"]  # Mother not in sample list

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        # Should still work but might not detect some patterns
        assert isinstance(patterns, list)

    def test_missing_genotypes(self, trio_pedigree):
        """Test handling of missing genotypes."""
        variant_row = {"father": "./.", "mother": "0/1", "child": "0/1"}
        sample_list = ["father", "mother", "child"]

        patterns = deduce_patterns_for_variant(variant_row, trio_pedigree, sample_list)
        # Should handle missing genotypes gracefully
        assert isinstance(patterns, list)

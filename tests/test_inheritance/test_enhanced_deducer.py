"""Test cases for enhanced pattern deduction logic."""

import pytest

from variantcentrifuge.inheritance.deducer import (
    check_dominant_pattern,
    check_mitochondrial_pattern,
    check_recessive_pattern,
    check_x_linked_patterns,
    deduce_patterns_for_variant,
    deduce_single_sample_patterns,
)

# Mark all tests in this module
pytestmark = [pytest.mark.inheritance, pytest.mark.deducer, pytest.mark.unit]


class TestSingleSampleDeduction:
    """Test single sample pattern deduction."""

    def test_single_sample_homozygous(self):
        """Test homozygous variant in single sample."""
        variant_row = {"sample1": "1/1"}
        sample_list = ["sample1"]

        patterns = deduce_single_sample_patterns(variant_row, sample_list)
        assert "homozygous" in patterns

    def test_single_sample_heterozygous(self):
        """Test heterozygous variant in single sample."""
        variant_row = {"sample1": "0/1"}
        sample_list = ["sample1"]

        patterns = deduce_single_sample_patterns(variant_row, sample_list)
        assert "unknown" in patterns

    def test_single_sample_reference(self):
        """Test reference genotype in single sample."""
        variant_row = {"sample1": "0/0"}
        sample_list = ["sample1"]

        patterns = deduce_single_sample_patterns(variant_row, sample_list)
        assert "reference" in patterns

    def test_single_sample_missing(self):
        """Test missing genotype in single sample."""
        variant_row = {"sample1": "./."}
        sample_list = ["sample1"]

        patterns = deduce_single_sample_patterns(variant_row, sample_list)
        assert patterns == ["unknown"]


class TestDeNovoDetection:
    """Test de novo pattern detection."""

    def test_classic_de_novo(self):
        """Test classic de novo: child has variant, parents don't."""
        variant_row = {"child": "0/1", "father": "0/0", "mother": "0/0"}
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

        patterns = deduce_patterns_for_variant(variant_row, pedigree_data, sample_list)
        assert "de_novo" in patterns

    def test_de_novo_candidate_missing_parent(self):
        """Test de novo candidate when parent genotype is missing."""
        variant_row = {"child": "0/1", "father": "./.", "mother": "0/0"}
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

        patterns = deduce_patterns_for_variant(variant_row, pedigree_data, sample_list)
        assert "de_novo_candidate" in patterns


class TestDominantPatterns:
    """Test autosomal dominant pattern detection."""

    def test_classic_dominant_from_father(self):
        """Test dominant inheritance from affected father."""
        variant_row = {"child": "0/1", "father": "0/1", "mother": "0/0"}
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
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = ["child", "father", "mother"]

        pattern = check_dominant_pattern("child", variant_row, pedigree_data, sample_list)
        assert pattern == "autosomal_dominant"

    def test_dominant_possible_incomplete_data(self):
        """Test possible dominant when parent data incomplete."""
        variant_row = {"child": "0/1", "father": "./.", "mother": "0/0"}
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

        pattern = check_dominant_pattern("child", variant_row, pedigree_data, sample_list)
        assert pattern == "autosomal_dominant_possible"


class TestRecessivePatterns:
    """Test autosomal recessive pattern detection."""

    def test_classic_recessive(self):
        """Test classic recessive: affected child homozygous, parents carriers."""
        variant_row = {"child": "1/1", "father": "0/1", "mother": "0/1"}
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

        pattern = check_recessive_pattern("child", variant_row, pedigree_data, sample_list)
        assert pattern == "autosomal_recessive"

    def test_recessive_possible_missing_parent(self):
        """Test possible recessive when parent data missing."""
        variant_row = {"child": "1/1", "father": "./.", "mother": "0/1"}
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

        pattern = check_recessive_pattern("child", variant_row, pedigree_data, sample_list)
        assert pattern == "autosomal_recessive_possible"


class TestXLinkedPatterns:
    """Test X-linked pattern detection."""

    def test_x_linked_recessive_male(self):
        """Test X-linked recessive in affected male."""
        variant_row = {
            "CHROM": "X",
            "son": "0/1",  # Hemizygous
            "father": "0/0",
            "mother": "0/1",  # Carrier
        }
        pedigree_data = {
            "son": {
                "sample_id": "son",
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
        sample_list = ["son", "father", "mother"]

        patterns = check_x_linked_patterns("son", variant_row, pedigree_data, sample_list)
        assert "x_linked_recessive" in patterns

    def test_x_linked_recessive_female(self):
        """Test X-linked recessive in affected female (homozygous)."""
        variant_row = {
            "CHROM": "X",
            "daughter": "1/1",
            "father": "0/1",  # Affected father
            "mother": "0/1",  # Carrier mother
        }
        pedigree_data = {
            "daughter": {
                "sample_id": "daughter",
                "father_id": "father",
                "mother_id": "mother",
                "sex": "2",
                "affected_status": "2",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
        }
        sample_list = ["daughter", "father", "mother"]

        patterns = check_x_linked_patterns("daughter", variant_row, pedigree_data, sample_list)
        assert "x_linked_recessive" in patterns

    def test_x_linked_violation_father_to_son(self):
        """Test that father-to-son transmission violates X-linked inheritance."""
        variant_row = {
            "CHROM": "X",
            "son": "0/1",
            "father": "0/1",  # Father has variant - violation!
            "mother": "0/0",
        }
        pedigree_data = {
            "son": {
                "sample_id": "son",
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
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            },
        }
        sample_list = ["son", "father", "mother"]

        patterns = check_x_linked_patterns("son", variant_row, pedigree_data, sample_list)
        assert "x_linked_recessive" not in patterns


class TestMitochondrialPatterns:
    """Test mitochondrial inheritance pattern detection."""

    def test_mitochondrial_maternal_transmission(self):
        """Test mitochondrial inheritance from mother."""
        variant_row = {"CHROM": "MT", "child": "1/1", "father": "0/0", "mother": "1/1"}
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
                "affected_status": "2",
            },
        }
        sample_list = ["child", "father", "mother"]

        result = check_mitochondrial_pattern("child", variant_row, pedigree_data, sample_list)
        assert result is True

    def test_mitochondrial_no_paternal_transmission(self):
        """Test that paternal transmission doesn't occur for mitochondrial."""
        variant_row = {
            "CHROM": "MT",
            "child": "0/0",
            "father": "1/1",  # Father has variant
            "mother": "0/0",  # Mother doesn't
        }
        pedigree_data = {
            "child": {
                "sample_id": "child",
                "father_id": "father",
                "mother_id": "mother",
                "affected_status": "1",
            },
            "father": {
                "sample_id": "father",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "2",
            },
            "mother": {
                "sample_id": "mother",
                "father_id": "0",
                "mother_id": "0",
                "affected_status": "1",
            },
        }
        sample_list = ["child", "father", "mother"]

        # Child doesn't have variant, which is consistent with no paternal transmission
        result = check_mitochondrial_pattern("child", variant_row, pedigree_data, sample_list)
        assert result is False  # Child doesn't have variant

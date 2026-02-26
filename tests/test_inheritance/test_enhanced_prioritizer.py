"""Test cases for enhanced pattern prioritizer with segregation support."""

from variantcentrifuge.inheritance.prioritizer import (
    PATTERN_PRIORITY,
    calculate_confidence,
    get_pattern_description,
    prioritize_patterns,
)


class TestPatternPrioritization:
    """Test pattern prioritization logic."""

    def test_basic_prioritization(self):
        """Test basic pattern prioritization without segregation."""
        patterns = ["unknown", "autosomal_dominant", "de_novo"]

        best_pattern, confidence = prioritize_patterns(patterns)

        # de_novo should win (highest priority)
        assert best_pattern == "de_novo"
        assert confidence > 0

    def test_prioritization_with_segregation(self):
        """Test prioritization with segregation data."""
        patterns = ["autosomal_dominant", "autosomal_recessive"]
        segregation_results = {
            "autosomal_dominant": (False, 0.2),  # Doesn't segregate well
            "autosomal_recessive": (True, 0.9),  # Segregates well
        }

        best_pattern, confidence = prioritize_patterns(
            patterns, segregation_results=segregation_results
        )

        # Recessive should win due to better segregation
        assert best_pattern == "autosomal_recessive"
        assert confidence > 0.5

    def test_single_sample_patterns(self):
        """Test prioritization of single sample patterns."""
        patterns = ["unknown", "homozygous"]

        best_pattern, _confidence = prioritize_patterns(patterns)

        # Homozygous has higher priority than unknown
        assert best_pattern == "homozygous"

    def test_compound_het_priority(self):
        """Test compound het patterns have high priority."""
        patterns = ["unknown", "compound_heterozygous_possible_no_pedigree"]

        best_pattern, _confidence = prioritize_patterns(patterns)

        # Compound het should win over unknown
        assert best_pattern == "compound_heterozygous_possible_no_pedigree"
        assert (
            PATTERN_PRIORITY["compound_heterozygous_possible_no_pedigree"]
            > PATTERN_PRIORITY["unknown"]
        )


class TestConfidenceCalculation:
    """Test confidence score calculation."""

    def test_single_pattern_confidence(self):
        """Test confidence when only one pattern."""
        pattern_scores = {"de_novo": 80}
        confidence = calculate_confidence(pattern_scores, "de_novo")

        assert 0 <= confidence <= 1
        assert confidence == 0.8  # 80/100

    def test_multiple_pattern_confidence(self):
        """Test confidence with multiple patterns."""
        pattern_scores = {"de_novo": 100, "autosomal_dominant": 60, "unknown": 10}
        confidence = calculate_confidence(pattern_scores, "de_novo")

        assert 0 <= confidence <= 1
        # Should have high confidence due to large separation
        assert confidence > 0.5

    def test_close_scores_low_confidence(self):
        """Test low confidence when scores are close."""
        pattern_scores = {"pattern1": 51, "pattern2": 50, "pattern3": 49}
        confidence = calculate_confidence(pattern_scores, "pattern1")

        assert confidence < 0.3  # Low confidence due to close scores


class TestPatternDescriptions:
    """Test pattern description generation."""

    def test_get_pattern_description(self):
        """Test getting human-readable descriptions."""
        assert "new mutation" in get_pattern_description("de_novo").lower()
        desc = get_pattern_description("compound_heterozygous").lower()
        assert "two different mutations" in desc or "compound" in desc
        assert "single sample" in get_pattern_description("homozygous").lower()
        assert get_pattern_description("nonexistent") == "Unknown inheritance pattern"


class TestPriorityValues:
    """Test that priority values are set correctly."""

    def test_pattern_priorities(self):
        """Test that pattern priorities follow expected order."""
        # De novo should have highest priority
        assert PATTERN_PRIORITY["de_novo"] > PATTERN_PRIORITY["compound_heterozygous"]
        assert PATTERN_PRIORITY["compound_heterozygous"] > PATTERN_PRIORITY["autosomal_recessive"]

        # Confirmed patterns should have higher priority than possible
        assert PATTERN_PRIORITY["de_novo"] > PATTERN_PRIORITY["de_novo_candidate"]
        assert (
            PATTERN_PRIORITY["compound_heterozygous"]
            > PATTERN_PRIORITY["compound_heterozygous_possible"]
        )

        # Unknown should have low priority
        assert PATTERN_PRIORITY["unknown"] < PATTERN_PRIORITY["homozygous"]

        # Reference and none should have lowest
        assert PATTERN_PRIORITY["reference"] < PATTERN_PRIORITY["unknown"]
        assert PATTERN_PRIORITY["none"] == 0

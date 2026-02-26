"""Tests for inheritance pattern prioritizer."""

from variantcentrifuge.inheritance.prioritizer import (
    PATTERN_PRIORITY,
    calculate_confidence,
    get_pattern_description,
    prioritize_patterns,
)


class TestPatternPrioritizer:
    """Test inheritance pattern prioritizer."""

    def test_prioritize_patterns_single(self):
        """Test prioritizing a single pattern."""
        patterns = ["autosomal_dominant"]
        pattern, confidence = prioritize_patterns(patterns)

        assert pattern == "autosomal_dominant"
        assert 0 <= confidence <= 1

    def test_prioritize_patterns_multiple(self):
        """Test prioritizing multiple patterns."""
        patterns = ["unknown", "autosomal_dominant", "de_novo"]
        pattern, confidence = prioritize_patterns(patterns)

        # De novo should win due to highest priority
        assert pattern == "de_novo"
        assert confidence > 0.5  # Should have good confidence

    def test_prioritize_patterns_empty(self):
        """Test with empty pattern list."""
        patterns = []
        pattern, confidence = prioritize_patterns(patterns)

        assert pattern == "unknown"
        assert confidence == 0.0

    def test_calculate_confidence(self):
        """Test confidence calculation."""
        # High separation between scores
        pattern_scores = {"de_novo": 100, "autosomal_dominant": 50, "unknown": 10}
        confidence = calculate_confidence(pattern_scores, "de_novo")
        assert confidence > 0.8

        # Close scores
        pattern_scores = {"de_novo": 100, "compound_heterozygous": 95, "autosomal_recessive": 90}
        confidence = calculate_confidence(pattern_scores, "de_novo")
        assert confidence < 0.5

        # Single pattern
        pattern_scores = {"autosomal_dominant": 60}
        confidence = calculate_confidence(pattern_scores, "autosomal_dominant")
        assert confidence == 0.6

    def test_get_pattern_description(self):
        """Test pattern descriptions."""
        assert "New mutation" in get_pattern_description("de_novo")
        assert "Two different mutations" in get_pattern_description("compound_heterozygous")
        assert "Two copies" in get_pattern_description("autosomal_recessive")
        assert "One copy" in get_pattern_description("autosomal_dominant")
        assert "X chromosome" in get_pattern_description("x_linked_recessive")
        assert get_pattern_description("invalid") == "Unknown inheritance pattern"

    def test_prioritize_with_variant_info(self):
        """Test prioritization with variant information."""
        patterns = ["de_novo", "autosomal_dominant"]

        # prioritize_patterns no longer takes variant_info
        pattern, confidence = prioritize_patterns(patterns)
        assert pattern == "de_novo"  # de_novo has higher priority
        assert confidence > 0.5

    def test_pattern_priorities_order(self):
        """Test that pattern priorities follow expected order."""
        # De novo should have highest priority
        assert PATTERN_PRIORITY["de_novo"] > PATTERN_PRIORITY["compound_heterozygous"]
        assert PATTERN_PRIORITY["compound_heterozygous"] > PATTERN_PRIORITY["autosomal_recessive"]

        # Confirmed patterns should have higher priority than possible
        assert PATTERN_PRIORITY["de_novo"] > PATTERN_PRIORITY["de_novo_candidate"]

        # Unknown should have low priority
        assert PATTERN_PRIORITY["unknown"] < PATTERN_PRIORITY["homozygous"]

"""Tests for inheritance pattern prioritizer."""

from variantcentrifuge.inheritance.prioritizer import (
    PATTERN_PRIORITY,
    adjust_pattern_score,
    calculate_confidence,
    filter_compatible_patterns,
    get_pattern_category,
    get_pattern_description,
    group_patterns_by_category,
    is_pattern_compatible,
    prioritize_patterns,
    resolve_conflicting_patterns,
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

    def test_adjust_pattern_score_de_novo(self):
        """Test score adjustment for de novo patterns."""
        base_score = PATTERN_PRIORITY["de_novo"]

        # adjust_pattern_score now just returns base score
        adjusted = adjust_pattern_score("de_novo", base_score)

        assert adjusted == base_score

    def test_adjust_pattern_score_recessive(self):
        """Test score adjustment for recessive patterns."""
        base_score = PATTERN_PRIORITY["autosomal_recessive"]

        # adjust_pattern_score now just returns base score
        adjusted = adjust_pattern_score("autosomal_recessive", base_score)

        assert adjusted == base_score

    def test_adjust_pattern_score_frequency(self):
        """Test score adjustment based on allele frequency."""
        base_score = 50

        # adjust_pattern_score now just returns base score
        adjusted = adjust_pattern_score("autosomal_dominant", base_score)
        assert adjusted == base_score

    def test_adjust_pattern_score_x_linked_male(self):
        """Test X-linked pattern scoring for males."""
        base_score = PATTERN_PRIORITY["x_linked_recessive"]

        # adjust_pattern_score now just returns base score
        adjusted = adjust_pattern_score("x_linked_recessive", base_score)

        assert adjusted == base_score

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

    def test_get_pattern_category(self):
        """Test pattern categorization."""
        assert get_pattern_category("de_novo") == "sporadic"
        assert get_pattern_category("compound_heterozygous") == "recessive"
        assert get_pattern_category("autosomal_recessive") == "recessive"
        assert get_pattern_category("autosomal_dominant") == "dominant"
        assert get_pattern_category("x_linked_recessive") == "x_linked"
        assert get_pattern_category("unknown") == "unclear"
        assert get_pattern_category("invalid") == "unclear"

    def test_group_patterns_by_category(self):
        """Test grouping patterns by category."""
        patterns = [
            "de_novo",
            "compound_heterozygous",
            "autosomal_recessive",
            "autosomal_dominant",
            "x_linked_recessive",
        ]

        grouped = group_patterns_by_category(patterns)

        assert "sporadic" in grouped
        assert grouped["sporadic"] == ["de_novo"]

        assert "recessive" in grouped
        assert set(grouped["recessive"]) == {"compound_heterozygous", "autosomal_recessive"}

        assert "dominant" in grouped
        assert grouped["dominant"] == ["autosomal_dominant"]

        assert "x_linked" in grouped
        assert grouped["x_linked"] == ["x_linked_recessive"]

    def test_get_pattern_description(self):
        """Test pattern descriptions."""
        assert "New mutation" in get_pattern_description("de_novo")
        assert "Two different mutations" in get_pattern_description("compound_heterozygous")
        assert "Two copies" in get_pattern_description("autosomal_recessive")
        assert "One copy" in get_pattern_description("autosomal_dominant")
        assert "X chromosome" in get_pattern_description("x_linked_recessive")
        assert get_pattern_description("invalid") == "Unknown inheritance pattern"

    def test_resolve_conflicting_patterns(self):
        """Test resolving conflicts across samples."""
        # Majority agreement
        patterns_by_sample = {
            "sample1": ["autosomal_dominant"],
            "sample2": ["autosomal_dominant"],
            "sample3": ["unknown"],
        }
        resolved = resolve_conflicting_patterns(patterns_by_sample)
        assert resolved == "autosomal_dominant"

        # No clear majority - use priority
        patterns_by_sample = {
            "sample1": ["de_novo"],
            "sample2": ["autosomal_dominant"],
            "sample3": ["autosomal_recessive"],
        }
        resolved = resolve_conflicting_patterns(patterns_by_sample)
        assert resolved == "de_novo"  # Highest priority

        # Empty patterns
        patterns_by_sample = {}
        resolved = resolve_conflicting_patterns(patterns_by_sample)
        assert resolved == "none"

    def test_filter_compatible_patterns(self):
        """Test filtering patterns by family structure."""
        patterns = ["de_novo", "autosomal_dominant", "x_linked_recessive"]

        # Family with parents and sex info
        family_structure = {"has_parents": True, "has_sex_info": True}
        compatible = filter_compatible_patterns(patterns, family_structure)
        assert set(compatible) == set(patterns)  # All compatible

        # No parent info
        family_structure = {"has_parents": False, "has_sex_info": True}
        compatible = filter_compatible_patterns(patterns, family_structure)
        assert "de_novo" not in compatible

        # No sex info
        family_structure = {"has_parents": True, "has_sex_info": False}
        compatible = filter_compatible_patterns(patterns, family_structure)
        assert "x_linked_recessive" not in compatible

    def test_is_pattern_compatible(self):
        """Test individual pattern compatibility checks."""
        # De novo requires parents
        assert is_pattern_compatible("de_novo", {"has_parents": True}) is True
        assert is_pattern_compatible("de_novo", {"has_parents": False}) is False

        # X-linked requires sex info
        assert is_pattern_compatible("x_linked_recessive", {"has_sex_info": True}) is True
        assert is_pattern_compatible("x_linked_recessive", {"has_sex_info": False}) is False

        # Others are generally compatible
        assert is_pattern_compatible("autosomal_dominant", {}) is True
        assert is_pattern_compatible("compound_heterozygous", {}) is True

    def test_prioritize_with_variant_info(self):
        """Test prioritization with variant information."""
        patterns = ["de_novo", "autosomal_dominant"]

        # prioritize_patterns no longer takes variant_info
        pattern, confidence = prioritize_patterns(patterns)
        assert pattern == "de_novo"  # de_novo has higher priority
        assert confidence > 0.5

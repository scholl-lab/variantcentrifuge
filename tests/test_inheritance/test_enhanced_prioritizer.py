"""Test cases for enhanced pattern prioritizer with segregation support."""

from variantcentrifuge.inheritance.prioritizer import (
    prioritize_patterns,
    adjust_pattern_score,
    calculate_confidence,
    get_pattern_category,
    group_patterns_by_category,
    get_pattern_description,
    resolve_conflicting_patterns,
    filter_compatible_patterns,
    PATTERN_PRIORITY,
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

        best_pattern, confidence = prioritize_patterns(patterns)

        # Homozygous has higher priority than unknown
        assert best_pattern == "homozygous"

    def test_compound_het_priority(self):
        """Test compound het patterns have high priority."""
        patterns = ["unknown", "compound_heterozygous_possible_no_pedigree"]

        best_pattern, confidence = prioritize_patterns(patterns)

        # Compound het should win over unknown
        assert best_pattern == "compound_heterozygous_possible_no_pedigree"
        assert (
            PATTERN_PRIORITY["compound_heterozygous_possible_no_pedigree"]
            > PATTERN_PRIORITY["unknown"]
        )


class TestScoreAdjustment:
    """Test pattern score adjustments."""

    def test_adjust_score_with_segregation(self):
        """Test score adjustment based on segregation."""
        base_score = 50
        
        # adjust_pattern_score now just returns base score
        score1 = adjust_pattern_score("pattern1", base_score)
        assert score1 == base_score
        
        score2 = adjust_pattern_score("pattern2", base_score)
        assert score2 == base_score

    def test_adjust_score_with_variant_info(self):
        """Test score adjustment based on variant information."""
        base_score = 50

        # adjust_pattern_score now just returns base score
        score = adjust_pattern_score("de_novo", base_score)
        assert score == base_score

        score = adjust_pattern_score("autosomal_recessive", base_score)
        assert score == base_score

        score = adjust_pattern_score("de_novo", base_score)
        assert score == base_score

    def test_adjust_score_with_sample_info(self):
        """Test score adjustment based on sample information."""
        base_score = 50

        # adjust_pattern_score now just returns base score
        score = adjust_pattern_score("x_linked_recessive", base_score)
        assert score == base_score

        score = adjust_pattern_score("autosomal_dominant", base_score)
        assert score == base_score


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


class TestPatternCategories:
    """Test pattern categorization."""

    def test_get_pattern_category(self):
        """Test getting category for patterns."""
        assert get_pattern_category("de_novo") == "sporadic"
        assert get_pattern_category("compound_heterozygous") == "recessive"
        assert get_pattern_category("autosomal_dominant") == "dominant"
        assert get_pattern_category("unknown") == "unclear"
        assert get_pattern_category("nonexistent_pattern") == "unclear"  # Default

    def test_group_patterns_by_category(self):
        """Test grouping patterns by category."""
        patterns = [
            "de_novo",
            "compound_heterozygous",
            "autosomal_recessive",
            "autosomal_dominant",
            "unknown",
        ]

        grouped = group_patterns_by_category(patterns)

        assert "sporadic" in grouped
        assert "de_novo" in grouped["sporadic"]
        assert "recessive" in grouped
        assert len(grouped["recessive"]) == 2  # compound_het and recessive
        assert "dominant" in grouped
        assert "unclear" in grouped


class TestPatternDescriptions:
    """Test pattern description generation."""

    def test_get_pattern_description(self):
        """Test getting human-readable descriptions."""
        assert "new mutation" in get_pattern_description("de_novo").lower()
        desc = get_pattern_description("compound_heterozygous").lower()
        assert "two different mutations" in desc or "compound" in desc
        assert "single sample" in get_pattern_description("homozygous").lower()
        assert get_pattern_description("nonexistent") == "Unknown inheritance pattern"


class TestConflictResolution:
    """Test resolving conflicting patterns across samples."""

    def test_resolve_unanimous_patterns(self):
        """Test when all samples agree on pattern."""
        patterns_by_sample = {
            "sample1": ["de_novo"],
            "sample2": ["de_novo"],
            "sample3": ["de_novo"],
        }

        result = resolve_conflicting_patterns(patterns_by_sample)
        assert result == "de_novo"

    def test_resolve_majority_patterns(self):
        """Test when majority agrees on pattern."""
        patterns_by_sample = {
            "sample1": ["autosomal_recessive"],
            "sample2": ["autosomal_recessive"],
            "sample3": ["unknown"],
            "sample4": ["autosomal_recessive"],
        }

        result = resolve_conflicting_patterns(patterns_by_sample)
        assert result == "autosomal_recessive"  # 3 out of 4

    def test_resolve_no_consensus(self):
        """Test when no clear consensus."""
        patterns_by_sample = {
            "sample1": ["de_novo"],
            "sample2": ["autosomal_dominant"],
            "sample3": ["autosomal_recessive"],
            "sample4": ["unknown"],
        }

        result = resolve_conflicting_patterns(patterns_by_sample)
        # Should use priority to resolve
        assert result in ["de_novo", "autosomal_dominant", "autosomal_recessive"]


class TestPatternFiltering:
    """Test filtering patterns by compatibility."""

    def test_filter_compatible_patterns(self):
        """Test filtering patterns based on family structure."""
        patterns = ["de_novo", "autosomal_dominant", "x_linked_recessive"]

        # Family with parents
        family_structure = {"has_parents": True, "has_sex_info": True}
        filtered = filter_compatible_patterns(patterns, family_structure)
        assert "de_novo" in filtered  # Requires parents
        assert "x_linked_recessive" in filtered  # Requires sex info

        # Single sample (no parents)
        family_structure = {"has_parents": False, "has_sex_info": False}
        filtered = filter_compatible_patterns(patterns, family_structure)
        assert "de_novo" not in filtered  # Can't determine de novo without parents
        assert "autosomal_dominant" in filtered  # Still compatible


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

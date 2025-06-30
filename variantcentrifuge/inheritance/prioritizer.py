"""
Pattern prioritizer for inheritance analysis.

This module prioritizes among multiple possible inheritance patterns
based on clinical significance and pattern reliability.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)

# Pattern priority scores (higher is better)
PATTERN_PRIORITY = {
    "de_novo": 100,  # Highest priority - often most clinically relevant
    "compound_heterozygous": 90,  # High priority - clear recessive mechanism
    "autosomal_recessive": 80,  # High priority - clear inheritance
    "x_linked_recessive": 75,  # High priority for males
    "x_linked_dominant": 70,  # Less common but significant
    "autosomal_dominant": 60,  # Common but clear pattern
    "mitochondrial": 50,  # Special inheritance
    "unknown": 10,  # Lowest priority
    "none": 0,  # No inheritance pattern
}

# Pattern categories for grouping
PATTERN_CATEGORIES = {
    "de_novo": "sporadic",
    "compound_heterozygous": "recessive",
    "autosomal_recessive": "recessive",
    "x_linked_recessive": "x_linked",
    "x_linked_dominant": "x_linked",
    "autosomal_dominant": "dominant",
    "mitochondrial": "maternal",
    "unknown": "unclear",
}


def prioritize_patterns(
    patterns: List[str],
    variant_info: Optional[Dict[str, Any]] = None,
    sample_info: Optional[Dict[str, Any]] = None,
) -> Tuple[str, float]:
    """
    Prioritize inheritance patterns and return the highest priority pattern.

    Args:
        patterns: List of possible inheritance patterns
        variant_info: Optional variant-specific information for scoring
        sample_info: Optional sample-specific information for scoring

    Returns:
        Tuple of (highest_priority_pattern, confidence_score)
    """
    if not patterns:
        return ("none", 0.0)

    # Calculate scores for each pattern
    pattern_scores = {}
    for pattern in patterns:
        base_score = PATTERN_PRIORITY.get(pattern, 0)

        # Adjust score based on additional information
        adjusted_score = adjust_pattern_score(pattern, base_score, variant_info, sample_info)

        pattern_scores[pattern] = adjusted_score

    # Get highest scoring pattern
    best_pattern = max(pattern_scores.items(), key=lambda x: x[1])

    # Calculate confidence based on score distribution
    confidence = calculate_confidence(pattern_scores, best_pattern[0])

    return (best_pattern[0], confidence)


def adjust_pattern_score(
    pattern: str,
    base_score: float,
    variant_info: Optional[Dict[str, Any]] = None,
    sample_info: Optional[Dict[str, Any]] = None,
) -> float:
    """
    Adjust pattern score based on additional evidence.

    Args:
        pattern: The inheritance pattern
        base_score: Base priority score
        variant_info: Variant-specific information
        sample_info: Sample-specific information

    Returns:
        Adjusted score
    """
    score = base_score

    if variant_info:
        # Boost de novo if variant is rare
        if pattern == "de_novo":
            if variant_info.get("is_rare", False):
                score += 20
            if variant_info.get("is_deleterious", False):
                score += 15

        # Boost recessive patterns for loss-of-function variants
        if pattern in ["autosomal_recessive", "compound_heterozygous"]:
            if variant_info.get("is_lof", False):
                score += 10

        # Consider allele frequency
        af = variant_info.get("allele_frequency", 1.0)
        if af < 0.001:  # Very rare
            score += 5
        elif af > 0.05:  # Common
            score -= 10

    if sample_info:
        # Boost X-linked patterns for males
        if pattern == "x_linked_recessive" and sample_info.get("sex") == "1":
            score += 10

        # Consider family history
        if sample_info.get("family_history", False):
            if pattern in ["autosomal_dominant", "autosomal_recessive"]:
                score += 5

    return max(score, 0)  # Don't allow negative scores


def calculate_confidence(pattern_scores: Dict[str, float], best_pattern: str) -> float:
    """
    Calculate confidence score based on pattern score distribution.

    Args:
        pattern_scores: Dictionary of pattern scores
        best_pattern: The highest scoring pattern

    Returns:
        Confidence score between 0 and 1
    """
    if not pattern_scores:
        return 0.0

    best_score = pattern_scores[best_pattern]

    # If only one pattern, confidence depends on its priority
    if len(pattern_scores) == 1:
        return min(best_score / 100.0, 1.0)

    # Calculate score separation
    other_scores = [score for pattern, score in pattern_scores.items() if pattern != best_pattern]

    if not other_scores:
        return 1.0

    second_best = max(other_scores)
    score_separation = best_score - second_best

    # Confidence based on separation and absolute score
    confidence = min((score_separation / 50.0) * (best_score / 100.0), 1.0)

    return round(confidence, 2)


def get_pattern_category(pattern: str) -> str:
    """
    Get the category for an inheritance pattern.

    Args:
        pattern: The inheritance pattern

    Returns:
        Pattern category
    """
    return PATTERN_CATEGORIES.get(pattern, "unclear")


def group_patterns_by_category(patterns: List[str]) -> Dict[str, List[str]]:
    """
    Group patterns by their categories.

    Args:
        patterns: List of inheritance patterns

    Returns:
        Dictionary mapping categories to patterns
    """
    grouped = {}
    for pattern in patterns:
        category = get_pattern_category(pattern)
        if category not in grouped:
            grouped[category] = []
        grouped[category].append(pattern)

    return grouped


def get_pattern_description(pattern: str) -> str:
    """
    Get a human-readable description of an inheritance pattern.

    Args:
        pattern: The inheritance pattern

    Returns:
        Description string
    """
    descriptions = {
        "de_novo": "New mutation not inherited from parents",
        "compound_heterozygous": "Two different mutations in the same gene, one from each parent",
        "autosomal_recessive": "Two copies of the mutation needed, one from each parent",
        "autosomal_dominant": "One copy of the mutation is sufficient to cause disease",
        "x_linked_recessive": "Mutation on X chromosome, primarily affects males",
        "x_linked_dominant": "Mutation on X chromosome affects both sexes",
        "mitochondrial": "Inherited through maternal lineage only",
        "unknown": "Inheritance pattern could not be determined",
        "none": "No inheritance pattern detected",
    }

    return descriptions.get(pattern, "Unknown inheritance pattern")


def resolve_conflicting_patterns(
    patterns_by_sample: Dict[str, List[str]], variant_info: Optional[Dict[str, Any]] = None
) -> str:
    """
    Resolve conflicting inheritance patterns across samples.

    Args:
        patterns_by_sample: Dictionary mapping sample IDs to their patterns
        variant_info: Optional variant information

    Returns:
        The most likely overall pattern
    """
    # Collect all patterns
    all_patterns = []
    for patterns in patterns_by_sample.values():
        all_patterns.extend(patterns)

    if not all_patterns:
        return "none"

    # Count pattern occurrences
    pattern_counts = {}
    for pattern in all_patterns:
        pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1

    # If one pattern dominates, use it
    total_count = len(all_patterns)
    for pattern, count in pattern_counts.items():
        if count / total_count > 0.5:  # More than 50% agreement
            return pattern

    # Otherwise, prioritize based on pattern priority
    unique_patterns = list(set(all_patterns))
    best_pattern, _ = prioritize_patterns(unique_patterns, variant_info)

    return best_pattern


def filter_compatible_patterns(patterns: List[str], family_structure: Dict[str, Any]) -> List[str]:
    """
    Filter patterns based on family structure compatibility.

    Args:
        patterns: List of potential patterns
        family_structure: Information about family structure

    Returns:
        List of patterns compatible with family structure
    """
    compatible = []

    for pattern in patterns:
        if is_pattern_compatible(pattern, family_structure):
            compatible.append(pattern)

    # If no patterns are compatible, return original list
    return compatible if compatible else patterns


def is_pattern_compatible(pattern: str, family_structure: Dict[str, Any]) -> bool:
    """
    Check if a pattern is compatible with family structure.

    Args:
        pattern: The inheritance pattern
        family_structure: Family structure information

    Returns:
        True if compatible
    """
    # De novo requires parents
    if pattern == "de_novo":
        return family_structure.get("has_parents", False)

    # X-linked patterns require known sex
    if pattern in ["x_linked_recessive", "x_linked_dominant"]:
        return family_structure.get("has_sex_info", False)

    # Compound het benefits from parent data
    if pattern == "compound_heterozygous":
        return True  # Can be inferred even without parents

    return True

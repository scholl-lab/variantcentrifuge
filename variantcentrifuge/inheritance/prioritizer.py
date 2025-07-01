"""
Pattern prioritizer for inheritance analysis.

This module prioritizes among multiple possible inheritance patterns
based on clinical significance and pattern reliability.
"""

import logging
from typing import List, Dict, Any, Optional, Tuple

logger = logging.getLogger(__name__)

# Pattern priority scores (higher is better) - matching variant-linker order
PATTERN_PRIORITY = {
    # Confirmed patterns with highest confidence
    "de_novo": 100,  # Highest priority - new mutations
    "compound_heterozygous": 90,  # Confirmed trans configuration
    "autosomal_recessive": 80,  # Both alleles affected
    "x_linked_recessive": 75,  # X-linked inheritance confirmed
    "x_linked_dominant": 70,  # X-linked dominant confirmed
    "autosomal_dominant": 60,  # Dominant inheritance confirmed
    "mitochondrial": 55,  # Maternal inheritance
    # Possible/candidate patterns
    "de_novo_candidate": 50,  # Possible de novo
    "compound_heterozygous_possible": 45,  # Possible compound het
    "compound_heterozygous_possible_no_pedigree": 44,  # Possible with no family data
    "compound_heterozygous_possible_missing_parents": 43,  # Missing parent data
    "compound_heterozygous_possible_missing_parent_genotypes": 42,  # Missing genotypes
    "autosomal_recessive_possible": 40,  # Possible recessive
    "x_linked_recessive_possible": 38,  # Possible X-linked recessive
    "x_linked_dominant_possible": 36,  # Possible X-linked dominant
    "autosomal_dominant_possible": 35,  # Possible dominant
    # Single sample patterns
    "homozygous": 25,  # Homozygous variant
    # Low priority patterns
    "non_mendelian": 12,  # Doesn't follow Mendelian rules
    "unknown": 10,  # Pattern unclear
    "reference": 5,  # No variant
    "none": 0,  # No pattern
}

# Pattern categories for grouping
PATTERN_CATEGORIES = {
    # Confirmed patterns
    "de_novo": "sporadic",
    "compound_heterozygous": "recessive",
    "autosomal_recessive": "recessive",
    "x_linked_recessive": "x_linked",
    "x_linked_dominant": "x_linked",
    "autosomal_dominant": "dominant",
    "mitochondrial": "maternal",
    # Possible patterns
    "de_novo_candidate": "sporadic",
    "compound_heterozygous_possible": "recessive",
    "compound_heterozygous_possible_no_pedigree": "recessive",
    "compound_heterozygous_possible_missing_parents": "recessive",
    "compound_heterozygous_possible_missing_parent_genotypes": "recessive",
    "autosomal_recessive_possible": "recessive",
    "x_linked_recessive_possible": "x_linked",
    "x_linked_dominant_possible": "x_linked",
    "autosomal_dominant_possible": "dominant",
    # Single sample patterns
    "homozygous": "recessive",
    "reference": "none",
    # Other patterns
    "non_mendelian": "complex",
    "unknown": "unclear",
    "none": "none",
}


def prioritize_patterns(
    patterns: List[str],
    variant_info: Optional[Dict[str, Any]] = None,
    sample_info: Optional[Dict[str, Any]] = None,
    segregation_results: Optional[Dict[str, Tuple[bool, float]]] = None,
) -> Tuple[str, float]:
    """
    Prioritize inheritance patterns and return the highest priority pattern.

    Parameters
    ----------
    patterns : List[str]
        List of possible inheritance patterns
    variant_info : Optional[Dict[str, Any]]
        Optional variant-specific information for scoring
    sample_info : Optional[Dict[str, Any]]
        Optional sample-specific information for scoring
    segregation_results : Optional[Dict[str, Tuple[bool, float]]]
        Optional segregation check results

    Returns
    -------
    Tuple[str, float]
        Tuple of (highest_priority_pattern, confidence_score)
    """
    if not patterns:
        return ("unknown", 0.0)

    # Filter patterns by segregation if results available
    if segregation_results:
        # Separate patterns that segregate vs those that don't
        segregating_patterns = []
        non_segregating_patterns = []

        for pattern in patterns:
            if pattern in segregation_results:
                segregates, seg_confidence = segregation_results[pattern]
                if segregates and seg_confidence > 0.5:
                    segregating_patterns.append(pattern)
                else:
                    non_segregating_patterns.append(pattern)
            else:
                # If no segregation data, include in segregating list
                segregating_patterns.append(pattern)

        # Strongly prefer segregating patterns
        if segregating_patterns:
            patterns = segregating_patterns
        elif non_segregating_patterns:
            # If no patterns segregate well, use all but penalize scores
            patterns = non_segregating_patterns

    # Calculate scores for each pattern
    pattern_scores = {}
    for pattern in patterns:
        base_score = PATTERN_PRIORITY.get(pattern, 0)

        # Adjust score based on additional information
        adjusted_score = adjust_pattern_score(
            pattern, base_score, variant_info, sample_info, segregation_results
        )

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
    segregation_results: Optional[Dict[str, Tuple[bool, float]]] = None,
) -> float:
    """
    Adjust pattern score based on additional evidence.

    Parameters
    ----------
    pattern : str
        The inheritance pattern
    base_score : float
        Base priority score
    variant_info : Optional[Dict[str, Any]]
        Variant-specific information
    sample_info : Optional[Dict[str, Any]]
        Sample-specific information
    segregation_results : Optional[Dict[str, Tuple[bool, float]]]
        Segregation check results

    Returns
    -------
    float
        Adjusted score
    """
    score = base_score

    # Apply segregation penalty/bonus
    if segregation_results and pattern in segregation_results:
        segregates, seg_confidence = segregation_results[pattern]
        if segregates:
            # Bonus for good segregation
            score += seg_confidence * 20
        else:
            # Penalty for poor segregation
            score -= (1 - seg_confidence) * 30

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

    Parameters
    ----------
    pattern_scores : Dict[str, float]
        Dictionary of pattern scores
    best_pattern : str
        The highest scoring pattern

    Returns
    -------
    float
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

    Parameters
    ----------
    pattern : str
        The inheritance pattern

    Returns
    -------
    str
        Pattern category
    """
    return PATTERN_CATEGORIES.get(pattern, "unclear")


def group_patterns_by_category(patterns: List[str]) -> Dict[str, List[str]]:
    """
    Group patterns by their categories.

    Parameters
    ----------
    patterns : List[str]
        List of inheritance patterns

    Returns
    -------
    Dict[str, List[str]]
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

    Parameters
    ----------
    pattern : str
        The inheritance pattern

    Returns
    -------
    str
        Description string
    """
    descriptions = {
        # Confirmed patterns
        "de_novo": "New mutation not inherited from parents",
        "compound_heterozygous": "Two different mutations in the same gene, one from each parent",
        "autosomal_recessive": "Two copies of the mutation needed, one from each parent",
        "autosomal_dominant": "One copy of the mutation is sufficient to cause disease",
        "x_linked_recessive": "Mutation on X chromosome, primarily affects males",
        "x_linked_dominant": "Mutation on X chromosome affects both sexes",
        "mitochondrial": "Inherited through maternal lineage only",
        # Possible/candidate patterns
        "de_novo_candidate": "Possible new mutation (parental data incomplete)",
        "compound_heterozygous_possible": "Possible compound heterozygous (phase uncertain)",
        "compound_heterozygous_possible_no_pedigree": (
            "Possible compound heterozygous (no family data)"
        ),
        "compound_heterozygous_possible_missing_parents": (
            "Possible compound heterozygous (missing parent data)"
        ),
        "compound_heterozygous_possible_missing_parent_genotypes": (
            "Possible compound heterozygous (missing genotypes)"
        ),
        "autosomal_recessive_possible": "Possible recessive inheritance (incomplete data)",
        "x_linked_recessive_possible": "Possible X-linked recessive (incomplete data)",
        "x_linked_dominant_possible": "Possible X-linked dominant (incomplete data)",
        "autosomal_dominant_possible": "Possible dominant inheritance (incomplete data)",
        # Single sample patterns
        "homozygous": "Homozygous variant (single sample)",
        "reference": "No variant present",
        # Other patterns
        "non_mendelian": "Does not follow Mendelian inheritance patterns",
        "unknown": "Inheritance pattern could not be determined",
        "none": "No inheritance pattern detected",
    }

    return descriptions.get(pattern, "Unknown inheritance pattern")


def resolve_conflicting_patterns(
    patterns_by_sample: Dict[str, List[str]], variant_info: Optional[Dict[str, Any]] = None
) -> str:
    """
    Resolve conflicting inheritance patterns across samples.

    Parameters
    ----------
    patterns_by_sample : Dict[str, List[str]]
        Dictionary mapping sample IDs to their patterns
    variant_info : Optional[Dict[str, Any]]
        Optional variant information

    Returns
    -------
    str
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

    Parameters
    ----------
    patterns : List[str]
        List of potential patterns
    family_structure : Dict[str, Any]
        Information about family structure

    Returns
    -------
    List[str]
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

    Parameters
    ----------
    pattern : str
        The inheritance pattern
    family_structure : Dict[str, Any]
        Family structure information

    Returns
    -------
    bool
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

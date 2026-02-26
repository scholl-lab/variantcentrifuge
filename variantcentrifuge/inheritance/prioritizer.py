"""
Pattern prioritizer for inheritance analysis.

This module prioritizes among multiple possible inheritance patterns
based on clinical significance and pattern reliability.
"""

import logging
from collections.abc import Mapping

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


def prioritize_patterns(
    patterns: list[str],
    segregation_results: dict[str, tuple[bool, float]] | None = None,
) -> tuple[str, float]:
    """
    Prioritize inheritance patterns and return the highest priority pattern.

    Parameters
    ----------
    patterns : List[str]
        List of possible inheritance patterns
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

        # Use segregating patterns if any exist, otherwise fall back to non-segregating
        if segregating_patterns:
            patterns = segregating_patterns
        elif non_segregating_patterns:
            patterns = non_segregating_patterns

    # Calculate scores for each pattern
    pattern_scores = {}
    for pattern in patterns:
        base_score = PATTERN_PRIORITY.get(pattern, 0)
        pattern_scores[pattern] = base_score

    # Get highest scoring pattern
    best_pattern = max(pattern_scores.items(), key=lambda x: x[1])

    # Calculate confidence based on score distribution
    confidence = calculate_confidence(pattern_scores, best_pattern[0])

    return (best_pattern[0], confidence)


def calculate_confidence(pattern_scores: Mapping[str, int | float], best_pattern: str) -> float:
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

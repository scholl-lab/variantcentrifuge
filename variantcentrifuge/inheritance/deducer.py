"""
Pattern deducer for Mendelian inheritance analysis.

This module contains the core logic for deducing inheritance patterns
based on genotypes in a family.
"""

import logging
from typing import Dict, List, Any
from ..genotype_utils import (
    is_variant,
    is_hom_alt,
    is_missing,
    could_be_de_novo,
)
from ..ped_reader import get_parents, is_affected

logger = logging.getLogger(__name__)


def deduce_patterns_for_variant(
    variant_row: Dict[str, Any], pedigree_data: Dict[str, Dict[str, Any]], sample_list: List[str]
) -> List[str]:
    """
    Deduce all possible inheritance patterns for a variant.

    Args:
        variant_row: Dictionary containing variant data with sample genotypes
        pedigree_data: Pedigree information dictionary
        sample_list: List of sample IDs to analyze

    Returns:
        List of possible inheritance patterns
    """
    patterns = []

    # Analyze each sample
    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        sample_patterns = deduce_patterns_for_sample(
            sample_id, variant_row, pedigree_data, sample_list
        )
        patterns.extend(sample_patterns)

    # Remove duplicates while preserving order
    seen = set()
    unique_patterns = []
    for pattern in patterns:
        if pattern not in seen:
            seen.add(pattern)
            unique_patterns.append(pattern)

    return unique_patterns


def deduce_patterns_for_sample(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> List[str]:
    """
    Deduce inheritance patterns for a specific sample.

    Args:
        sample_id: The sample to analyze
        variant_row: Variant data with genotypes
        pedigree_data: Pedigree information
        sample_list: All available samples

    Returns:
        List of possible patterns for this sample
    """
    patterns = []

    # Get sample genotype
    sample_gt = variant_row.get(sample_id, "./.")

    # Skip if no variant
    if not is_variant(sample_gt):
        return patterns

    # Get parent information
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Check for de novo
    if father_id and mother_id and father_id in sample_list and mother_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        mother_gt = variant_row.get(mother_id, "./.")

        if could_be_de_novo(sample_gt, father_gt, mother_gt):
            patterns.append("de_novo")

    # Check for dominant patterns
    if is_dominant_pattern(sample_id, variant_row, pedigree_data, sample_list):
        patterns.append("autosomal_dominant")

    # Check for recessive patterns
    if is_recessive_pattern(sample_id, variant_row, pedigree_data, sample_list):
        patterns.append("autosomal_recessive")

    # Check for X-linked patterns
    sex = pedigree_data[sample_id].get("sex", "0")
    if sex in ["1", "2"]:  # Known sex
        if is_x_linked_pattern(sample_id, variant_row, pedigree_data, sample_list):
            if sex == "1":  # Male
                patterns.append("x_linked_recessive")
            else:  # Female
                patterns.append("x_linked_dominant")

    # If no specific pattern found but variant present
    if not patterns and is_variant(sample_gt):
        patterns.append("unknown")

    return patterns


def is_dominant_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> bool:
    """
    Check if variant follows autosomal dominant pattern.

    Criteria:
    - Affected individuals have the variant
    - At least one parent has the variant (unless de novo)
    - Transmitted from affected parent to affected offspring
    """
    sample_gt = variant_row.get(sample_id, "./.")

    # Must have variant
    if not is_variant(sample_gt):
        return False

    # Check if affected
    if not is_affected(sample_id, pedigree_data):
        return False

    # Check parents
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    if father_id and father_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        if is_variant(father_gt) and is_affected(father_id, pedigree_data):
            return True

    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        if is_variant(mother_gt) and is_affected(mother_id, pedigree_data):
            return True

    return False


def is_recessive_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> bool:
    """
    Check if variant follows autosomal recessive pattern.

    Criteria:
    - Affected individuals are homozygous for the variant
    - Parents are typically carriers (heterozygous)
    - Both parents must carry at least one copy
    """
    sample_gt = variant_row.get(sample_id, "./.")

    # Must be homozygous alt
    if not is_hom_alt(sample_gt):
        return False

    # Should be affected
    if not is_affected(sample_id, pedigree_data):
        return False

    # Check parents
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Both parents should be carriers or affected
    father_carrier = False
    mother_carrier = False

    if father_id and father_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        father_carrier = is_variant(father_gt)

    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        mother_carrier = is_variant(mother_gt)

    # Both parents must carry the variant
    return father_carrier and mother_carrier


def is_x_linked_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> bool:
    """
    Check if variant follows X-linked pattern.

    Criteria:
    - Males: Hemizygous (one copy is sufficient)
    - Females: Need two copies for recessive, one for dominant
    - Typically transmitted from carrier mothers to affected sons
    """
    sample_gt = variant_row.get(sample_id, "./.")
    sex = pedigree_data[sample_id].get("sex", "0")

    # Must have variant
    if not is_variant(sample_gt):
        return False

    # Should be affected
    if not is_affected(sample_id, pedigree_data):
        return False

    # Male-specific logic
    if sex == "1":  # Male
        # Males need only one copy (hemizygous)
        # Check if mother is carrier
        father_id, mother_id = get_parents(sample_id, pedigree_data)
        if mother_id and mother_id in sample_list:
            mother_gt = variant_row.get(mother_id, "./.")
            return is_variant(mother_gt)

    # Female-specific logic
    elif sex == "2":  # Female
        # For X-linked recessive, females need to be homozygous
        # For X-linked dominant, heterozygous is sufficient
        return is_variant(sample_gt)

    return False


def check_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
    pattern: str,
) -> bool:
    """
    Check if variant segregates according to the specified pattern.

    Args:
        variant_row: Variant data with genotypes
        pedigree_data: Pedigree information
        sample_list: All available samples
        pattern: The inheritance pattern to check

    Returns:
        True if segregation is consistent with pattern
    """
    # Count affected and unaffected individuals with variants
    affected_with_variant = 0
    affected_without_variant = 0
    unaffected_with_variant = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        gt = variant_row.get(sample_id, "./.")
        if is_missing(gt):
            continue

        has_variant = is_variant(gt)
        affected = is_affected(sample_id, pedigree_data)

        if affected and has_variant:
            affected_with_variant += 1
        elif affected and not has_variant:
            affected_without_variant += 1
        elif not affected and has_variant:
            unaffected_with_variant += 1

    # Pattern-specific checks
    if pattern == "autosomal_dominant":
        # Most affected should have variant, few unaffected should have it
        return affected_with_variant > 0 and affected_without_variant == 0

    elif pattern == "autosomal_recessive":
        # Affected should be homozygous, unaffected can be carriers
        return affected_with_variant > 0 and affected_without_variant == 0

    return True


def get_inheritance_info(
    sample_id: str,
    patterns: List[str],
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Get detailed inheritance information for a sample.

    Args:
        sample_id: The sample ID
        patterns: List of deduced patterns
        variant_row: Variant data
        pedigree_data: Pedigree information

    Returns:
        Dictionary with detailed inheritance information
    """
    info = {
        "sample_id": sample_id,
        "patterns": patterns,
        "genotype": variant_row.get(sample_id, "./."),
        "affected_status": is_affected(sample_id, pedigree_data),
    }

    # Add parent information
    father_id, mother_id = get_parents(sample_id, pedigree_data)
    if father_id:
        info["father_genotype"] = variant_row.get(father_id, "./.")
        info["father_affected"] = is_affected(father_id, pedigree_data)
    if mother_id:
        info["mother_genotype"] = variant_row.get(mother_id, "./.")
        info["mother_affected"] = is_affected(mother_id, pedigree_data)

    return info

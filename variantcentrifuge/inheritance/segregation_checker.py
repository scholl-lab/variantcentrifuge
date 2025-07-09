"""
Segregation checker for inheritance analysis.

This module validates whether inheritance patterns segregate correctly
with the affected status across family members.
"""

import logging
from typing import Any, Dict, List, Tuple

from ..genotype_utils import is_het, is_hom_alt, is_missing, is_variant
from ..ped_reader import get_parents, is_affected

logger = logging.getLogger(__name__)


def check_pattern_segregation(
    pattern: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """
    Check if a variant segregates according to a specific inheritance pattern.

    Parameters
    ----------
    pattern : str
        The inheritance pattern to check
    variant_row : Dict[str, Any]
        Variant data with genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        Available samples

    Returns
    -------
    Tuple[bool, float]
        Tuple of (segregates_correctly, confidence_score)
    """
    if pattern == "de_novo":
        return check_de_novo_segregation(variant_row, pedigree_data, sample_list)
    elif pattern == "autosomal_dominant":
        return check_dominant_segregation(variant_row, pedigree_data, sample_list)
    elif pattern == "autosomal_recessive":
        return check_recessive_segregation(variant_row, pedigree_data, sample_list)
    elif pattern == "compound_heterozygous":
        return check_compound_het_segregation(variant_row, pedigree_data, sample_list)
    elif pattern in ["x_linked_recessive", "x_linked_dominant"]:
        return check_x_linked_segregation(pattern, variant_row, pedigree_data, sample_list)
    elif pattern == "mitochondrial":
        return check_mitochondrial_segregation(variant_row, pedigree_data, sample_list)
    else:
        # Unknown patterns default to True with low confidence
        return (True, 0.3)


def check_de_novo_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """Check de novo segregation: variant in child but not in parents."""
    violations = 0
    checks = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        father_id, mother_id = get_parents(sample_id, pedigree_data)

        # Only check samples with both parents available
        if father_id and mother_id and father_id in sample_list and mother_id in sample_list:
            sample_gt = variant_row.get(sample_id, "./.")
            father_gt = variant_row.get(father_id, "./.")
            mother_gt = variant_row.get(mother_id, "./.")

            if (
                not is_missing(sample_gt)
                and not is_missing(father_gt)
                and not is_missing(mother_gt)
            ):
                checks += 1

                # De novo: child has variant, parents don't
                if is_variant(sample_gt):
                    if is_variant(father_gt) or is_variant(mother_gt):
                        violations += 1

    if checks == 0:
        return (True, 0.5)  # No data to check

    confidence = 1.0 - (violations / checks)
    segregates = violations == 0

    return (segregates, confidence)


def check_dominant_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """
    Check autosomal dominant segregation.

    - All affected individuals should have the variant
    - Unaffected individuals should not have the variant (with some tolerance)
    """
    affected_with_variant = 0
    affected_without_variant = 0
    unaffected_with_variant = 0
    unaffected_without_variant = 0

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
        else:  # not affected and not has_variant
            unaffected_without_variant += 1

    total_affected = affected_with_variant + affected_without_variant
    total_unaffected = unaffected_with_variant + unaffected_without_variant

    if total_affected == 0:
        return (False, 0.0)  # No affected individuals

    # Perfect segregation: all affected have variant, no unaffected have it
    if affected_without_variant == 0 and unaffected_with_variant == 0:
        return (True, 1.0)

    # Calculate confidence based on violations
    affected_rate = affected_with_variant / total_affected if total_affected > 0 else 0
    unaffected_rate = unaffected_without_variant / total_unaffected if total_unaffected > 0 else 1

    confidence = (affected_rate + unaffected_rate) / 2

    # Segregates if most affected have it and few unaffected have it
    segregates = affected_rate > 0.8 and (unaffected_with_variant / (total_unaffected + 1)) < 0.2

    return (segregates, confidence)


def check_recessive_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """
    Check autosomal recessive segregation.

    - Affected individuals should be homozygous
    - Parents of affected should be carriers
    - Unaffected can be carriers or wild-type
    """
    violations = 0
    checks = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        gt = variant_row.get(sample_id, "./.")
        if is_missing(gt):
            continue

        affected = is_affected(sample_id, pedigree_data)

        if affected:
            checks += 1
            # Affected must be homozygous
            if not is_hom_alt(gt):
                violations += 1

            # Check parents if available
            father_id, mother_id = get_parents(sample_id, pedigree_data)
            if father_id and mother_id and father_id in sample_list and mother_id in sample_list:
                father_gt = variant_row.get(father_id, "./.")
                mother_gt = variant_row.get(mother_id, "./.")

                # Both parents should carry at least one variant allele
                if not is_missing(father_gt) and not is_missing(mother_gt):
                    checks += 1
                    if not (is_variant(father_gt) and is_variant(mother_gt)):
                        violations += 1

    if checks == 0:
        return (False, 0.0)

    confidence = 1.0 - (violations / checks)
    segregates = violations == 0

    return (segregates, confidence)


def check_compound_het_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """
    Check compound heterozygous segregation.

    Note: This requires checking multiple variants in the same gene,
    which is handled separately in comp_het.py
    """
    # For a single variant, we can only check if affected individuals
    # are heterozygous (necessary but not sufficient)
    affected_het_count = 0
    affected_total = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        if is_affected(sample_id, pedigree_data):
            gt = variant_row.get(sample_id, "./.")
            if not is_missing(gt):
                affected_total += 1
                if is_het(gt):
                    affected_het_count += 1

    if affected_total == 0:
        return (False, 0.0)

    # All affected should be heterozygous for compound het
    confidence = affected_het_count / affected_total
    segregates = affected_het_count == affected_total

    return (segregates, confidence)


def check_x_linked_segregation(
    pattern: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """Check X-linked segregation patterns."""
    violations = 0
    checks = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        gt = variant_row.get(sample_id, "./.")
        if is_missing(gt):
            continue

        sex = pedigree_data[sample_id].get("sex", "0")
        affected = is_affected(sample_id, pedigree_data)
        has_variant = is_variant(gt)

        checks += 1

        if pattern == "x_linked_recessive":
            if sex == "1":  # Male
                # Affected males must have variant (hemizygous)
                if affected and not has_variant:
                    violations += 1
                # Unaffected males should not have variant
                elif not affected and has_variant:
                    violations += 1
            elif sex == "2":  # Female
                # Affected females must be homozygous
                if affected and not is_hom_alt(gt):
                    violations += 1

        elif pattern == "x_linked_dominant":
            # Both males and females can be affected with one copy
            if affected and not has_variant:
                violations += 1

    if checks == 0:
        return (True, 0.5)

    confidence = 1.0 - (violations / checks)
    segregates = violations < (checks * 0.2)  # Allow 20% violations

    return (segregates, confidence)


def check_mitochondrial_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Tuple[bool, float]:
    """Check mitochondrial inheritance: maternal transmission only."""
    violations = 0
    checks = 0

    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        father_id, mother_id = get_parents(sample_id, pedigree_data)

        # Check maternal transmission
        if mother_id and mother_id in sample_list:
            sample_gt = variant_row.get(sample_id, "./.")
            mother_gt = variant_row.get(mother_id, "./.")

            if not is_missing(sample_gt) and not is_missing(mother_gt):
                checks += 1
                # If mother has variant, child should have it
                if is_variant(mother_gt) and not is_variant(sample_gt):
                    violations += 1

        # Check no paternal transmission
        if father_id and father_id in sample_list:
            sample_gt = variant_row.get(sample_id, "./.")
            father_gt = variant_row.get(father_id, "./.")

            if not is_missing(sample_gt) and not is_missing(father_gt):
                checks += 1
                # If father has variant, child should NOT have it (unless from mother)
                if is_variant(father_gt) and is_variant(sample_gt):
                    if not (
                        mother_id
                        and mother_id in sample_list
                        and is_variant(variant_row.get(mother_id, "./."))
                    ):
                        violations += 1

    if checks == 0:
        return (True, 0.5)

    confidence = 1.0 - (violations / checks)
    segregates = violations == 0

    return (segregates, confidence)


def calculate_segregation_score(
    patterns: List[str],
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Dict[str, Tuple[bool, float]]:
    """
    Calculate segregation scores for multiple patterns.

    Parameters
    ----------
    patterns : List[str]
        List of inheritance patterns to check
    variant_row : Dict[str, Any]
        Variant data with genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        Available samples

    Returns
    -------
    Dict[str, Tuple[bool, float]]
        Dictionary mapping patterns to (segregates, confidence) tuples
    """
    results = {}

    for pattern in patterns:
        segregates, confidence = check_pattern_segregation(
            pattern, variant_row, pedigree_data, sample_list
        )
        results[pattern] = (segregates, confidence)

    return results

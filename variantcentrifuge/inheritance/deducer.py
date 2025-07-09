"""
Pattern deducer for Mendelian inheritance analysis.

This module contains the core logic for deducing inheritance patterns
based on genotypes in a family.
"""

import logging
from typing import Any, Dict, List, Optional

from ..genotype_utils import could_be_de_novo, is_het, is_hom_alt, is_missing, is_ref, is_variant
from ..ped_reader import get_parents, is_affected

logger = logging.getLogger(__name__)


def deduce_patterns_for_variant(
    variant_row: Dict[str, Any], pedigree_data: Dict[str, Dict[str, Any]], sample_list: List[str]
) -> List[str]:
    """
    Deduce all possible inheritance patterns for a variant.

    Parameters
    ----------
    variant_row : Dict[str, Any]
        Dictionary containing variant data with sample genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information dictionary
    sample_list : List[str]
        List of sample IDs to analyze

    Returns
    -------
    List[str]
        List of possible inheritance patterns
    """
    patterns = []

    # Handle single sample case (no pedigree data)
    if not pedigree_data or len(sample_list) == 1:
        return deduce_single_sample_patterns(variant_row, sample_list)

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


def deduce_single_sample_patterns(variant_row: Dict[str, Any], sample_list: List[str]) -> List[str]:
    """
    Deduce patterns for single sample analysis (no family data).

    Parameters
    ----------
    variant_row : Dict[str, Any]
        Variant data with genotypes
    sample_list : List[str]
        List of sample IDs (typically just one)

    Returns
    -------
    List[str]
        List of basic patterns based on genotype
    """
    patterns = []

    for sample_id in sample_list:
        gt = variant_row.get(sample_id, "./.")

        if is_missing(gt):
            continue
        elif is_ref(gt):
            patterns.append("reference")
        elif is_hom_alt(gt):
            patterns.append("homozygous")
        elif is_het(gt):
            # For single sample, heterozygous pattern is unknown
            # Will be resolved by compound het analysis if multiple variants in same gene
            patterns.append("unknown")

    return patterns if patterns else ["unknown"]


def deduce_patterns_for_sample(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> List[str]:
    """
    Deduce inheritance patterns for a specific sample.

    Parameters
    ----------
    sample_id : str
        The sample to analyze
    variant_row : Dict[str, Any]
        Variant data with genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        All available samples

    Returns
    -------
    List[str]
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
        elif is_missing(father_gt) or is_missing(mother_gt):
            # If parent genotypes are missing, it could be de novo
            if is_variant(sample_gt) and is_affected(sample_id, pedigree_data):
                patterns.append("de_novo_candidate")

    # Check for dominant patterns
    dom_result = check_dominant_pattern(sample_id, variant_row, pedigree_data, sample_list)
    if dom_result:
        patterns.append(dom_result)

    # Check for recessive patterns
    rec_result = check_recessive_pattern(sample_id, variant_row, pedigree_data, sample_list)
    if rec_result:
        patterns.append(rec_result)

    # Check for X-linked patterns if on X chromosome
    chrom = str(variant_row.get("CHROM", "")).upper()
    if chrom in ["X", "CHRX", "23"]:
        x_patterns = check_x_linked_patterns(sample_id, variant_row, pedigree_data, sample_list)
        patterns.extend(x_patterns)

    # Check for mitochondrial if on MT chromosome
    if chrom in ["MT", "M", "CHRM", "CHRMT"]:
        if check_mitochondrial_pattern(sample_id, variant_row, pedigree_data, sample_list):
            patterns.append("mitochondrial")

    # If no specific pattern found but variant present
    if not patterns and is_variant(sample_gt):
        if is_affected(sample_id, pedigree_data):
            patterns.append("unknown")
        else:
            patterns.append("carrier")

    return patterns


def check_dominant_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Optional[str]:
    """
    Check if variant follows autosomal dominant pattern.

    Returns specific pattern type or None.
    """
    sample_gt = variant_row.get(sample_id, "./.")

    # Must have variant
    if not is_variant(sample_gt):
        return None

    # Check if affected
    if not is_affected(sample_id, pedigree_data):
        return None

    # Check parents
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    father_has_variant = False
    mother_has_variant = False
    parent_data_complete = True

    if father_id and father_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        if is_missing(father_gt):
            parent_data_complete = False
        elif is_variant(father_gt):
            father_has_variant = True

    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        if is_missing(mother_gt):
            parent_data_complete = False
        elif is_variant(mother_gt):
            mother_has_variant = True

    # Classic dominant: at least one parent affected and has variant
    if father_has_variant or mother_has_variant:
        if (father_has_variant and is_affected(father_id, pedigree_data)) or (
            mother_has_variant and is_affected(mother_id, pedigree_data)
        ):
            return "autosomal_dominant"
        else:
            # Parent has variant but not affected - incomplete penetrance?
            return "autosomal_dominant_possible"

    # If parent data incomplete
    if not parent_data_complete:
        return "autosomal_dominant_possible"

    return None


def check_recessive_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Optional[str]:
    """
    Check if variant follows autosomal recessive pattern.

    Returns specific pattern type or None.
    """
    sample_gt = variant_row.get(sample_id, "./.")

    # Must be homozygous alt for classic recessive
    if not is_hom_alt(sample_gt):
        return None

    # Should be affected
    if not is_affected(sample_id, pedigree_data):
        return None

    # Check parents
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Both parents should be carriers or affected
    father_carrier = None
    mother_carrier = None
    parent_data_complete = True

    if father_id and father_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        if is_missing(father_gt):
            parent_data_complete = False
        else:
            father_carrier = is_variant(father_gt)

    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        if is_missing(mother_gt):
            parent_data_complete = False
        else:
            mother_carrier = is_variant(mother_gt)

    # Classic recessive: both parents are carriers
    if father_carrier and mother_carrier:
        # Check if parents are affected
        father_affected = is_affected(father_id, pedigree_data) if father_id else False
        mother_affected = is_affected(mother_id, pedigree_data) if mother_id else False

        # If both parents are carriers but not affected, classic recessive
        if not father_affected and not mother_affected:
            return "autosomal_recessive"
        else:
            # One or both parents affected - still recessive but unusual
            return "autosomal_recessive"

    # If we have incomplete parent data but what we have is consistent
    if (
        not parent_data_complete
        and (father_carrier is None or father_carrier)
        and (mother_carrier is None or mother_carrier)
    ):
        return "autosomal_recessive_possible"

    return None


def check_x_linked_patterns(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> List[str]:
    """
    Check if variant follows X-linked patterns.

    Returns list of possible X-linked patterns.
    """
    patterns = []
    sample_gt = variant_row.get(sample_id, "./.")
    sex = pedigree_data[sample_id].get("sex", "0")
    affected = is_affected(sample_id, pedigree_data)

    # Must have variant
    if not is_variant(sample_gt):
        return patterns

    # Get parent information
    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Check X-linked recessive
    xlr_consistent = check_x_linked_recessive(sample_id, variant_row, pedigree_data, sample_list)
    if xlr_consistent:
        patterns.append(xlr_consistent)

    # Check X-linked dominant
    xld_consistent = check_x_linked_dominant(sample_id, variant_row, pedigree_data, sample_list)
    if xld_consistent:
        patterns.append(xld_consistent)

    # If affected but patterns unclear, add possible
    if affected and not patterns:
        if sex == "1":  # Male
            patterns.append("x_linked_recessive_possible")
        else:  # Female
            if is_hom_alt(sample_gt):
                patterns.append("x_linked_recessive_possible")
            else:
                patterns.append("x_linked_dominant_possible")

    return patterns


def check_x_linked_recessive(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Optional[str]:
    """Check for X-linked recessive pattern with sophisticated logic."""
    sample_gt = variant_row.get(sample_id, "./.")
    sex = pedigree_data[sample_id].get("sex", "0")
    affected = is_affected(sample_id, pedigree_data)

    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Male logic
    if sex == "1":  # Male
        if affected and is_variant(sample_gt):
            # Check mother is carrier
            if mother_id and mother_id in sample_list:
                mother_gt = variant_row.get(mother_id, "./.")
                if is_missing(mother_gt):
                    return "x_linked_recessive_possible"
                elif is_variant(mother_gt):
                    # Check father doesn't have it (males can't pass X to sons)
                    if father_id and father_id in sample_list:
                        father_gt = variant_row.get(father_id, "./.")
                        if is_variant(father_gt):
                            return None  # Violation: father can't pass X to son
                    return "x_linked_recessive"
                else:
                    return None  # Mother must be carrier
            else:
                return "x_linked_recessive_possible"

    # Female logic
    elif sex == "2":  # Female
        if affected and is_hom_alt(sample_gt):
            # Both parents must contribute
            father_has = None
            mother_has = None

            if father_id and father_id in sample_list:
                father_gt = variant_row.get(father_id, "./.")
                if not is_missing(father_gt):
                    father_has = is_variant(father_gt)

            if mother_id and mother_id in sample_list:
                mother_gt = variant_row.get(mother_id, "./.")
                if not is_missing(mother_gt):
                    mother_has = is_variant(mother_gt)

            if father_has and mother_has:
                return "x_linked_recessive"
            elif father_has is None or mother_has is None:
                return "x_linked_recessive_possible"

    return None


def check_x_linked_dominant(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> Optional[str]:
    """Check for X-linked dominant pattern."""
    sample_gt = variant_row.get(sample_id, "./.")
    sex = pedigree_data[sample_id].get("sex", "0")
    affected = is_affected(sample_id, pedigree_data)

    if not affected or not is_variant(sample_gt):
        return None

    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Check inheritance patterns
    if father_id and father_id in sample_list:
        father_gt = variant_row.get(father_id, "./.")
        father_affected = is_affected(father_id, pedigree_data)

        # If affected father has variant and child is female, consistent with XLD
        if father_affected and is_variant(father_gt) and sex == "2":
            return "x_linked_dominant"

    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        mother_affected = is_affected(mother_id, pedigree_data)

        # If affected mother has variant, consistent with XLD
        if mother_affected and is_variant(mother_gt):
            return "x_linked_dominant"

    # If we can't determine parent status clearly
    if affected and is_variant(sample_gt):
        return "x_linked_dominant_possible"

    return None


def check_mitochondrial_pattern(
    sample_id: str,
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
) -> bool:
    """Check for mitochondrial inheritance pattern."""
    sample_gt = variant_row.get(sample_id, "./.")

    if not is_variant(sample_gt):
        return False

    father_id, mother_id = get_parents(sample_id, pedigree_data)

    # Check maternal transmission
    if mother_id and mother_id in sample_list:
        mother_gt = variant_row.get(mother_id, "./.")
        if not is_missing(mother_gt):
            # If mother has it, consistent with mitochondrial
            return is_variant(mother_gt)

    # If no mother data, could still be mitochondrial
    return True


def check_segregation(
    variant_row: Dict[str, Any],
    pedigree_data: Dict[str, Dict[str, Any]],
    sample_list: List[str],
    pattern: str,
) -> bool:
    """
    Check if variant segregates according to the specified pattern.

    Parameters
    ----------
    variant_row : Dict[str, Any]
        Variant data with genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information
    sample_list : List[str]
        All available samples
    pattern : str
        The inheritance pattern to check

    Returns
    -------
    bool
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

    Parameters
    ----------
    sample_id : str
        The sample ID
    patterns : List[str]
        List of deduced patterns
    variant_row : Dict[str, Any]
        Variant data
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information

    Returns
    -------
    Dict[str, Any]
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

"""
Vectorized inheritance pattern deduction using NumPy boolean masks.

This module provides a vectorized implementation of Pass 1 (pattern deduction)
for inheritance analysis. Instead of using df.apply() which calls Python
per-row, it operates on entire columns using NumPy boolean masks for 10-100x
speedup.

Performance benefits:
- NumPy vectorized operations eliminate Python interpreter overhead
- Boolean masks operate on entire arrays at once
- Memory-efficient int8 encoding for genotypes
"""

import logging
import time
from typing import Any

import numpy as np
import pandas as pd

from ..ped_reader import get_parents, is_affected
from .comp_het_vectorized import encode_genotypes

logger = logging.getLogger(__name__)


def vectorized_deduce_patterns(
    df: pd.DataFrame, pedigree_data: dict[str, dict[str, Any]], sample_list: list[str]
) -> list[list[str]]:
    """
    Vectorized inheritance pattern deduction for ALL variants at once.

    This function replaces the df.apply(deduce_patterns_for_variant) approach
    with NumPy boolean mask operations for significant performance improvement.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing variant data with sample genotypes
    pedigree_data : Dict[str, Dict[str, Any]]
        Pedigree information dictionary
    sample_list : List[str]
        List of sample IDs to analyze

    Returns
    -------
    List[List[str]]
        List of pattern lists (one per variant), matching the output
        format of the original deduce_patterns_for_variant
    """
    t_start = time.monotonic()
    n_variants = len(df)

    # Initialize result structure - list of lists matching original format
    patterns_per_variant: list[list[str]] = [[] for _ in range(n_variants)]

    # Handle single sample / no pedigree case
    if not pedigree_data or len(sample_list) == 1:
        return _deduce_single_sample_patterns_vectorized(df, sample_list, n_variants)

    # Step 1: Build genotype matrix (n_variants x n_samples)
    # This encodes ALL sample genotypes ONCE before pattern checking
    t_matrix_start = time.monotonic()
    gt_matrix, sample_to_idx = _build_genotype_matrix(df, sample_list)
    t_matrix = time.monotonic() - t_matrix_start
    logger.info(
        f"Vectorized deducer: genotype matrix built in {t_matrix:.3f}s "
        f"({n_variants} variants x {len(sample_list)} samples, "
        f"matrix shape={gt_matrix.shape}, dtype={gt_matrix.dtype})"
    )

    # Step 2: Per-sample vectorized pattern deduction
    for sample_id in sample_list:
        if sample_id not in pedigree_data:
            continue

        # Get sample index in matrix
        if sample_id not in sample_to_idx:
            continue

        sample_idx = sample_to_idx[sample_id]
        sample_gts = gt_matrix[:, sample_idx]

        # Create has_variant mask (gt > 0 means het or hom_alt)
        has_variant_mask = sample_gts > 0

        # Get parent information
        father_id, mother_id = get_parents(sample_id, pedigree_data)
        affected = is_affected(sample_id, pedigree_data)

        # Get parent genotypes if available
        father_gts = None
        mother_gts = None
        if father_id and father_id in sample_to_idx:
            father_gts = gt_matrix[:, sample_to_idx[father_id]]
        if mother_id and mother_id in sample_to_idx:
            mother_gts = gt_matrix[:, sample_to_idx[mother_id]]

        # Track which patterns were added by this sample (for fallback logic)
        patterns_before_sample = [len(p) for p in patterns_per_variant]

        # a. De novo check (vectorized)
        _check_de_novo_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            sample_list,
            has_variant_mask,
            affected,
        )

        # b. Dominant check (vectorized)
        _check_dominant_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            pedigree_data,
            has_variant_mask,
            affected,
        )

        # c. Recessive check (vectorized)
        _check_recessive_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            affected,
        )

        # d. X-linked check (vectorized)
        chrom_array = df["CHROM"].values
        # Handle categorical dtypes by converting to string
        if isinstance(chrom_array.dtype, pd.CategoricalDtype):
            chrom_array = chrom_array.astype(str)
        else:
            chrom_array = chrom_array.astype(str)

        _check_x_linked_vectorized(
            patterns_per_variant,
            sample_id,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            pedigree_data,
            sample_list,
            chrom_array,
            has_variant_mask,
            affected,
        )

        # e. Mitochondrial check (vectorized)
        _check_mitochondrial_vectorized(
            patterns_per_variant,
            sample_gts,
            mother_gts,
            mother_id,
            sample_list,
            chrom_array,
            has_variant_mask,
        )

        # f. Fallback for variants with no patterns FROM THIS SAMPLE
        # This matches original logic: if this sample added no patterns but has variant,
        # add "unknown" (affected) or "carrier" (unaffected)
        _apply_fallback_vectorized(
            patterns_per_variant, sample_gts, has_variant_mask, affected, patterns_before_sample
        )

    t_patterns = time.monotonic() - t_start - t_matrix
    logger.info(
        f"Vectorized deducer: pattern checks completed in {t_patterns:.3f}s "
        f"({len(sample_list)} samples x {n_variants} variants)"
    )

    # Step 3: Deduplication per variant (preserve order)
    for i in range(n_variants):
        if patterns_per_variant[i]:
            # Remove duplicates while preserving order
            seen = set()
            unique_patterns = []
            for pattern in patterns_per_variant[i]:
                if pattern not in seen:
                    seen.add(pattern)
                    unique_patterns.append(pattern)
            patterns_per_variant[i] = unique_patterns

    t_total = time.monotonic() - t_start
    n_with_patterns = sum(1 for p in patterns_per_variant if p)
    logger.info(
        f"Vectorized deducer complete in {t_total:.3f}s: "
        f"{n_with_patterns}/{n_variants} variants have patterns "
        f"(matrix={t_matrix:.3f}s, patterns={t_patterns:.3f}s)"
    )

    return patterns_per_variant


def _build_genotype_matrix(
    df: pd.DataFrame, sample_list: list[str]
) -> tuple[np.ndarray, dict[str, int]]:
    """
    Build a genotype matrix (n_variants x n_samples) with int8 encoding.

    Parameters
    ----------
    df : pd.DataFrame
        Variant DataFrame
    sample_list : List[str]
        List of sample IDs

    Returns
    -------
    Tuple[np.ndarray, Dict[str, int]]
        - Genotype matrix (n_variants x n_samples)
        - Mapping from sample_id to column index
    """
    n_variants = len(df)
    n_samples = len(sample_list)

    # Initialize matrix with -1 (missing)
    gt_matrix = np.full((n_variants, n_samples), -1, dtype=np.int8)

    # Build sample to index mapping
    sample_to_idx = {}

    for i, sample_id in enumerate(sample_list):
        sample_to_idx[sample_id] = i

        # Encode genotypes for this sample if column exists
        if sample_id in df.columns:
            gt_matrix[:, i] = encode_genotypes(df[sample_id])

    return gt_matrix, sample_to_idx


def _deduce_single_sample_patterns_vectorized(
    df: pd.DataFrame, sample_list: list[str], n_variants: int
) -> list[list[str]]:
    """
    Vectorized single sample pattern deduction (no pedigree data).

    Parameters
    ----------
    df : pd.DataFrame
        Variant DataFrame
    sample_list : List[str]
        List of sample IDs
    n_variants : int
        Number of variants

    Returns
    -------
    List[List[str]]
        Pattern list per variant
    """
    patterns_per_variant: list[list[str]] = [[] for _ in range(n_variants)]

    for sample_id in sample_list:
        if sample_id not in df.columns:
            continue

        # Encode genotypes
        sample_gts = encode_genotypes(df[sample_id])

        # Create masks for each genotype type
        ref_mask = sample_gts == 0
        hom_alt_mask = sample_gts == 2
        het_mask = sample_gts == 1

        # Apply patterns based on genotype
        # Reference: append "reference"
        ref_indices = np.where(ref_mask)[0]
        for idx in ref_indices:
            patterns_per_variant[idx].append("reference")

        # Homozygous alt: append "homozygous"
        hom_alt_indices = np.where(hom_alt_mask)[0]
        for idx in hom_alt_indices:
            patterns_per_variant[idx].append("homozygous")

        # Heterozygous: append "unknown"
        het_indices = np.where(het_mask)[0]
        for idx in het_indices:
            patterns_per_variant[idx].append("unknown")

    # For variants with no pattern from any sample, append "unknown"
    for i in range(n_variants):
        if not patterns_per_variant[i]:
            patterns_per_variant[i].append("unknown")

    return patterns_per_variant


def _check_de_novo_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    sample_list: list[str],
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """
    Vectorized de novo pattern detection.

    De novo: child has variant, both parents are ref (0/0).
    De novo candidate: child has variant, at least one parent has missing GT.
    """
    # Both parents must be in sample list
    if not (father_id and mother_id and father_id in sample_list and mother_id in sample_list):
        return

    if father_gts is None or mother_gts is None:
        return

    # Classic de novo: child has variant, both parents ref
    de_novo_mask = has_variant_mask & (father_gts == 0) & (mother_gts == 0)
    de_novo_indices = np.where(de_novo_mask)[0]
    for idx in de_novo_indices:
        patterns_per_variant[idx].append("de_novo")

    # De novo candidate: child has variant, at least one parent missing GT, child affected
    if affected:
        de_novo_candidate_mask = has_variant_mask & ((father_gts == -1) | (mother_gts == -1))
        de_novo_candidate_indices = np.where(de_novo_candidate_mask)[0]
        for idx in de_novo_candidate_indices:
            patterns_per_variant[idx].append("de_novo_candidate")


def _check_dominant_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    pedigree_data: dict[str, dict[str, Any]],
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """
    Vectorized autosomal dominant pattern detection.

    Classic AD: affected child has variant, at least one affected parent has variant.
    AD possible: affected child has variant, parent has variant but not affected,
                 or parent data incomplete.
    """
    # Must be affected
    if not affected:
        return

    # Must have variant
    if not has_variant_mask.any():
        return

    # Check father contribution
    father_has_variant_mask = np.zeros(len(sample_gts), dtype=bool)
    father_affected = np.bool_(False)
    father_data_complete = np.ones(len(sample_gts), dtype=bool)

    if father_id and father_gts is not None:
        father_affected = np.bool_(is_affected(father_id, pedigree_data))
        father_has_variant_mask = father_gts > 0
        father_data_complete = father_gts != -1

    # Check mother contribution
    mother_has_variant_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_affected = np.bool_(False)
    mother_data_complete = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_gts is not None:
        mother_affected = np.bool_(is_affected(mother_id, pedigree_data))
        mother_has_variant_mask = mother_gts > 0
        mother_data_complete = mother_gts != -1

    # Classic dominant: at least one parent affected and has variant
    classic_ad_mask = has_variant_mask & (
        (father_has_variant_mask & father_affected) | (mother_has_variant_mask & mother_affected)
    )
    classic_ad_indices = np.where(classic_ad_mask)[0]
    for idx in classic_ad_indices:
        patterns_per_variant[idx].append("autosomal_dominant")

    # Parent has variant but not affected - incomplete penetrance
    parent_variant_not_affected_mask = has_variant_mask & (
        (father_has_variant_mask & ~father_affected) | (mother_has_variant_mask & ~mother_affected)
    )
    # Exclude indices already classified as classic AD
    parent_variant_not_affected_mask[classic_ad_indices] = False
    parent_variant_not_affected_indices = np.where(parent_variant_not_affected_mask)[0]
    for idx in parent_variant_not_affected_indices:
        patterns_per_variant[idx].append("autosomal_dominant_possible")

    # Incomplete parent data
    incomplete_data_mask = has_variant_mask & (~father_data_complete | ~mother_data_complete)
    # Exclude indices already classified
    incomplete_data_mask[classic_ad_indices] = False
    incomplete_data_mask[parent_variant_not_affected_indices] = False
    incomplete_data_indices = np.where(incomplete_data_mask)[0]
    for idx in incomplete_data_indices:
        patterns_per_variant[idx].append("autosomal_dominant_possible")


def _check_recessive_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    affected: bool,
) -> None:
    """
    Vectorized autosomal recessive pattern detection.

    Classic AR: affected child is hom_alt (1/1), both parents are carriers (have variant).
    AR possible: affected child is hom_alt, incomplete parent data consistent with AR.
    """
    # Must be affected
    if not affected:
        return

    # Must be homozygous alt
    hom_alt_mask = sample_gts == 2

    if not hom_alt_mask.any():
        return

    # Check parent carrier status
    father_carrier_mask = np.zeros(len(sample_gts), dtype=bool)
    father_data_complete = np.ones(len(sample_gts), dtype=bool)

    if father_id and father_gts is not None:
        father_carrier_mask = father_gts > 0
        father_data_complete = father_gts != -1

    mother_carrier_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_data_complete = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_gts is not None:
        mother_carrier_mask = mother_gts > 0
        mother_data_complete = mother_gts != -1

    # Classic recessive: both parents are carriers
    classic_ar_mask = hom_alt_mask & father_carrier_mask & mother_carrier_mask
    classic_ar_indices = np.where(classic_ar_mask)[0]
    for idx in classic_ar_indices:
        patterns_per_variant[idx].append("autosomal_recessive")

    # Incomplete parent data consistent with AR
    # (parent data incomplete but what we have is consistent)
    incomplete_ar_mask = (
        hom_alt_mask
        & ~classic_ar_mask
        & (~father_data_complete | ~mother_data_complete)
        & (father_carrier_mask | ~father_data_complete)
        & (mother_carrier_mask | ~mother_data_complete)
    )
    incomplete_ar_indices = np.where(incomplete_ar_mask)[0]
    for idx in incomplete_ar_indices:
        patterns_per_variant[idx].append("autosomal_recessive_possible")


def _check_x_linked_vectorized(
    patterns_per_variant: list[list[str]],
    sample_id: str,
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    pedigree_data: dict[str, dict[str, Any]],
    sample_list: list[str],
    chrom_array: np.ndarray,
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """
    Vectorized X-linked pattern detection.

    X-linked recessive (males): affected male has variant, mother is carrier.
    X-linked recessive (females): affected female is hom_alt, both parents have variant.
    X-linked dominant: affected individual has variant, affected parent has variant.
    """
    # Create X chromosome mask
    x_mask = np.isin(chrom_array, ["X", "CHRX", "23", "x", "chrX", "chrx"])

    # Must be on X chromosome
    if not x_mask.any():
        return

    # Get sex from pedigree
    sex = pedigree_data[sample_id].get("sex", "0")

    # Male X-linked recessive logic
    if sex == "1":  # Male
        _check_xlr_male_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            sample_list,
            x_mask,
            has_variant_mask,
            affected,
        )

    # Female X-linked logic
    elif sex == "2":  # Female
        _check_xlr_female_vectorized(
            patterns_per_variant,
            sample_gts,
            father_gts,
            mother_gts,
            father_id,
            mother_id,
            sample_list,
            x_mask,
            affected,
        )

    # X-linked dominant (both sexes)
    _check_xld_vectorized(
        patterns_per_variant,
        sample_id,
        sample_gts,
        father_gts,
        mother_gts,
        father_id,
        mother_id,
        pedigree_data,
        sex,
        x_mask,
        has_variant_mask,
        affected,
    )

    # Fallback for affected but unclear patterns
    if affected:
        _check_x_fallback_vectorized(
            patterns_per_variant, sample_gts, sex, x_mask, has_variant_mask
        )


def _check_xlr_male_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    sample_list: list[str],
    x_mask: np.ndarray,
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """X-linked recessive for males: affected + has variant + mother carrier."""
    if not affected:
        return

    # Male must have variant on X
    male_xlr_base_mask = x_mask & has_variant_mask

    if not male_xlr_base_mask.any():
        return

    # Check mother is carrier
    if mother_id and mother_id in sample_list and mother_gts is not None:
        mother_missing_mask = mother_gts == -1
        mother_has_variant_mask = mother_gts > 0

        # Father check: males can't pass X to sons, so father must NOT have variant
        father_violates_mask = np.zeros(len(sample_gts), dtype=bool)
        if father_id and father_id in sample_list and father_gts is not None:
            father_violates_mask = father_gts > 0

        # Classic XLR: mother has variant, father doesn't (or no father data)
        classic_xlr_mask = male_xlr_base_mask & mother_has_variant_mask & ~father_violates_mask
        classic_xlr_indices = np.where(classic_xlr_mask)[0]
        for idx in classic_xlr_indices:
            patterns_per_variant[idx].append("x_linked_recessive")

        # XLR possible: mother missing genotype
        xlr_possible_mask = male_xlr_base_mask & mother_missing_mask & ~father_violates_mask
        # Exclude classic XLR indices
        xlr_possible_mask[classic_xlr_indices] = False
        xlr_possible_indices = np.where(xlr_possible_mask)[0]
        for idx in xlr_possible_indices:
            patterns_per_variant[idx].append("x_linked_recessive_possible")
    else:
        # No mother data - XLR possible
        xlr_possible_indices = np.where(male_xlr_base_mask)[0]
        for idx in xlr_possible_indices:
            patterns_per_variant[idx].append("x_linked_recessive_possible")


def _check_xlr_female_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    sample_list: list[str],
    x_mask: np.ndarray,
    affected: bool,
) -> None:
    """X-linked recessive for females: affected + hom_alt + both parents have variant."""
    if not affected:
        return

    # Female must be hom_alt on X
    female_xlr_base_mask = x_mask & (sample_gts == 2)

    if not female_xlr_base_mask.any():
        return

    # Both parents must contribute
    father_has_mask = np.zeros(len(sample_gts), dtype=bool)
    father_missing_mask = np.ones(len(sample_gts), dtype=bool)

    if father_id and father_id in sample_list and father_gts is not None:
        father_has_mask = father_gts > 0
        father_missing_mask = father_gts == -1

    mother_has_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_missing_mask = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_id in sample_list and mother_gts is not None:
        mother_has_mask = mother_gts > 0
        mother_missing_mask = mother_gts == -1

    # Classic XLR: both parents have variant
    classic_xlr_mask = female_xlr_base_mask & father_has_mask & mother_has_mask
    classic_xlr_indices = np.where(classic_xlr_mask)[0]
    for idx in classic_xlr_indices:
        patterns_per_variant[idx].append("x_linked_recessive")

    # XLR possible: at least one parent has missing data
    xlr_possible_mask = (
        female_xlr_base_mask & (father_missing_mask | mother_missing_mask) & ~classic_xlr_mask
    )
    xlr_possible_indices = np.where(xlr_possible_mask)[0]
    for idx in xlr_possible_indices:
        patterns_per_variant[idx].append("x_linked_recessive_possible")


def _check_xld_vectorized(
    patterns_per_variant: list[list[str]],
    sample_id: str,
    sample_gts: np.ndarray,
    father_gts: np.ndarray | None,
    mother_gts: np.ndarray | None,
    father_id: str | None,
    mother_id: str | None,
    pedigree_data: dict[str, dict[str, Any]],
    sex: str,
    x_mask: np.ndarray,
    has_variant_mask: np.ndarray,
    affected: bool,
) -> None:
    """X-linked dominant: affected individual has variant, affected parent has variant."""
    if not affected:
        return

    # Must have variant on X
    xld_base_mask = x_mask & has_variant_mask

    if not xld_base_mask.any():
        return

    # Check father (XLD: affected father with variant can pass to daughters)
    if father_id and father_gts is not None:
        father_affected = np.bool_(is_affected(father_id, pedigree_data))
        father_has_variant_mask = father_gts > 0

        # If affected father has variant and child is female, XLD
        if father_affected and sex == "2":
            xld_father_mask = xld_base_mask & father_has_variant_mask
            xld_father_indices = np.where(xld_father_mask)[0]
            for idx in xld_father_indices:
                patterns_per_variant[idx].append("x_linked_dominant")

    # Check mother (XLD: affected mother with variant can pass to any child)
    if mother_id and mother_gts is not None:
        mother_affected = np.bool_(is_affected(mother_id, pedigree_data))
        mother_has_variant_mask = mother_gts > 0

        if mother_affected:
            xld_mother_mask = xld_base_mask & mother_has_variant_mask
            xld_mother_indices = np.where(xld_mother_mask)[0]
            for idx in xld_mother_indices:
                patterns_per_variant[idx].append("x_linked_dominant")

    # XLD possible: affected with variant but can't determine parent status
    # Check if we already added XLD
    xld_possible_mask = xld_base_mask.copy()
    for idx in np.where(xld_base_mask)[0]:
        if "x_linked_dominant" in patterns_per_variant[idx]:
            xld_possible_mask[idx] = False

    xld_possible_indices = np.where(xld_possible_mask)[0]
    for idx in xld_possible_indices:
        patterns_per_variant[idx].append("x_linked_dominant_possible")


def _check_x_fallback_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    sex: str,
    x_mask: np.ndarray,
    has_variant_mask: np.ndarray,
) -> None:
    """
    Fallback for affected individuals on X with unclear patterns.

    This adds "possible" patterns when affected but no clear XLR/XLD detected.
    """
    # Only apply if no X-linked patterns already added
    x_variant_mask = x_mask & has_variant_mask

    for idx in np.where(x_variant_mask)[0]:
        # Check if any X-linked pattern already present
        existing_patterns = patterns_per_variant[idx]
        has_x_pattern = any(
            p
            for p in existing_patterns
            if p.startswith("x_linked_recessive") or p.startswith("x_linked_dominant")
        )

        if not has_x_pattern:
            # Add sex-specific fallback
            if sex == "1":  # Male
                patterns_per_variant[idx].append("x_linked_recessive_possible")
            else:  # Female
                if sample_gts[idx] == 2:  # hom_alt
                    patterns_per_variant[idx].append("x_linked_recessive_possible")
                else:
                    patterns_per_variant[idx].append("x_linked_dominant_possible")


def _check_mitochondrial_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    mother_gts: np.ndarray | None,
    mother_id: str | None,
    sample_list: list[str],
    chrom_array: np.ndarray,
    has_variant_mask: np.ndarray,
) -> None:
    """
    Vectorized mitochondrial pattern detection.

    Mitochondrial: maternal inheritance (mother has variant, or no mother data).
    """
    # Create mitochondrial chromosome mask
    mt_mask = np.isin(chrom_array, ["MT", "M", "CHRM", "CHRMT", "mt", "chrM", "chrMT"])

    # Must be on MT chromosome and have variant
    mt_variant_mask = mt_mask & has_variant_mask

    if not mt_variant_mask.any():
        return

    # Check maternal transmission
    if mother_id and mother_id in sample_list and mother_gts is not None:
        # Mother must have variant for mitochondrial
        mother_has_variant_mask = mother_gts > 0

        # Maternal variant confirms mitochondrial
        maternal_mt_mask = mt_variant_mask & mother_has_variant_mask
        maternal_mt_indices = np.where(maternal_mt_mask)[0]
        for idx in maternal_mt_indices:
            patterns_per_variant[idx].append("mitochondrial")

        # If mother doesn't have variant, not mitochondrial (skip these)
        # If mother GT missing, still could be mitochondrial (apply below)
    else:
        # No mother data - could still be mitochondrial
        mt_indices = np.where(mt_variant_mask)[0]
        for idx in mt_indices:
            patterns_per_variant[idx].append("mitochondrial")


def _apply_fallback_vectorized(
    patterns_per_variant: list[list[str]],
    sample_gts: np.ndarray,
    has_variant_mask: np.ndarray,
    affected: bool,
    patterns_before_sample: list[int],
) -> None:
    """
    Apply fallback patterns for variants where THIS SAMPLE has variant but added no patterns.

    This matches original per-sample logic: if a sample has a variant but contributed
    no specific patterns, add "unknown" (affected) or "carrier" (unaffected).

    Parameters
    ----------
    patterns_per_variant : List[List[str]]
        Pattern lists per variant
    sample_gts : np.ndarray
        Sample genotypes
    has_variant_mask : np.ndarray
        Boolean mask where sample has variant
    affected : bool
        Whether sample is affected
    patterns_before_sample : List[int]
        Number of patterns per variant before this sample's analysis
    """
    # Find variants where sample has variant
    variant_indices = np.where(has_variant_mask)[0]

    for idx in variant_indices:
        # Check if THIS SAMPLE added any patterns
        # (compare current count to count before this sample)
        patterns_added_by_sample = len(patterns_per_variant[idx]) - patterns_before_sample[idx]

        if patterns_added_by_sample == 0:
            # This sample has variant but added no patterns - apply fallback
            if affected:
                patterns_per_variant[idx].append("unknown")
            else:
                patterns_per_variant[idx].append("carrier")

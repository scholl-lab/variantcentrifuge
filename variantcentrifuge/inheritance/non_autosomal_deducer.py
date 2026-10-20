"""
Vectorized non-autosomal (X-linked and mitochondrial) inheritance pattern detection.

Cohesive module for non-autosomal vectorized inheritance pattern analysis.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..ped_reader import is_affected

# GRCh38 PAR1/PAR2 (1-based, inclusive); PAR variants are inherited like autosomes.
# ponytail: the genome build is not passed into inheritance analysis, so GRCh38 PARs are used
# for every build. On GRCh37, PAR2 stays X-linked (old behaviour) and chrX:2,699,521-2,781,479
# counts as PAR; thread links.normalize_reference(config["reference"]) here to make it exact.
_PAR_GRCH38 = ((10_001, 2_781_479), (155_701_383, 156_030_895))
# PLINK numbering: 23=X, 24=Y, 25=XY (PAR, autosome-like), 26=MT
_NON_AUTOSOMAL = ["X", "Y", "M", "MT", "23", "24", "26"]

__all__ = [
    "_check_mitochondrial_vectorized",
    "_check_x_fallback_vectorized",
    "_check_x_linked_vectorized",
    "_check_xld_vectorized",
    "_check_xlr_female_vectorized",
    "_check_xlr_male_vectorized",
    "autosomal_mask",
]


def autosomal_mask(chrom: Any, pos: Any) -> np.ndarray:
    """Return True for autosomal and X pseudoautosomal (PAR) variants."""
    chrom_norm = (
        pd.Series(chrom, dtype=str).str.strip().str.upper().str.replace("^CHR", "", regex=True)
    )
    pos_num = pd.to_numeric(pd.Series(pos), errors="coerce").to_numpy(dtype=float)
    in_par = np.zeros(len(pos_num), dtype=bool)
    for start, end in _PAR_GRCH38:
        in_par |= (pos_num >= start) & (pos_num <= end)
    is_x = chrom_norm.isin(["X", "23"]).to_numpy()
    return np.asarray(~chrom_norm.isin(_NON_AUTOSOMAL).to_numpy() | (is_x & in_par), dtype=bool)


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
    chrom_clean = np.char.upper(np.char.strip(np.asarray(chrom_array, dtype=str)))
    x_mask = np.isin(chrom_clean, ["X", "CHRX", "23", "CHR23"])

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
            father_violates_mask = father_gts > 0  # type: ignore[assignment]

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
        father_has_mask = father_gts > 0  # type: ignore[assignment]
        father_missing_mask = father_gts == -1  # type: ignore[assignment]

    mother_has_mask = np.zeros(len(sample_gts), dtype=bool)
    mother_missing_mask = np.ones(len(sample_gts), dtype=bool)

    if mother_id and mother_id in sample_list and mother_gts is not None:
        mother_has_mask = mother_gts > 0  # type: ignore[assignment]
        mother_missing_mask = mother_gts == -1  # type: ignore[assignment]

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

    xld_found = np.zeros(len(sample_gts), dtype=bool)

    # Check father (XLD: affected father with variant can pass to daughters)
    if father_id and father_gts is not None:
        father_affected = np.bool_(is_affected(father_id, pedigree_data))
        father_has_variant_mask = father_gts > 0

        # If affected father has variant and child is female, XLD
        if father_affected and sex == "2":
            xld_father_mask = xld_base_mask & father_has_variant_mask
            xld_found |= xld_father_mask
            xld_father_indices = np.where(xld_father_mask)[0]
            for idx in xld_father_indices:
                patterns_per_variant[idx].append("x_linked_dominant")

    # Check mother (XLD: affected mother with variant can pass to any child)
    if mother_id and mother_gts is not None:
        mother_affected = np.bool_(is_affected(mother_id, pedigree_data))
        mother_has_variant_mask = mother_gts > 0

        if mother_affected:
            xld_mother_mask = xld_base_mask & mother_has_variant_mask
            xld_found |= xld_mother_mask
            xld_mother_indices = np.where(xld_mother_mask)[0]
            for idx in xld_mother_indices:
                patterns_per_variant[idx].append("x_linked_dominant")

    # XLD possible: no definite XLD for THIS sample (other samples' labels must not count)
    xld_possible_mask = xld_base_mask & ~xld_found
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
    chrom_clean = np.char.upper(np.char.strip(np.asarray(chrom_array, dtype=str)))
    mt_mask = np.isin(chrom_clean, ["MT", "M", "CHRM", "CHRMT", "26", "CHR26"])

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
    else:
        # No mother data - could still be mitochondrial
        mt_indices = np.where(mt_variant_mask)[0]
        for idx in mt_indices:
            patterns_per_variant[idx].append("mitochondrial")

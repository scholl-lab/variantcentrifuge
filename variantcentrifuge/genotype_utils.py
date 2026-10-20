"""
Genotype utility functions for inheritance analysis.

This module provides functions to parse and analyze genotype strings
commonly found in VCF files.
"""


def parse_genotype(gt: str) -> tuple[int | None, ...]:
    """
    Parse a genotype string into allele indices.

    Parameters
    ----------
    gt : str
        Genotype string (e.g., "0/1", "1|1", "0", "1", "./.", ".")

    Returns
    -------
    tuple[int | None, ...]
        Tuple of allele indices where None indicates missing
    """
    if not gt or gt in ("./.", "."):
        return (None,) if gt == "." else (None, None)

    # Haploid call (single allele without / or |)
    if "/" not in gt and "|" not in gt:
        try:
            return (int(gt),)
        except ValueError:
            return (None, None)

    # Handle both / and | separators
    separator = "|" if "|" in gt else "/"
    parts = gt.split(separator)

    if len(parts) != 2:
        return (None, None)

    try:
        allele1 = None if parts[0] == "." else int(parts[0])
        allele2 = None if parts[1] == "." else int(parts[1])
        return (allele1, allele2)
    except (ValueError, IndexError):
        return (None, None)


def is_het(gt: str) -> bool:
    """Check if genotype is heterozygous (0/1 or 1/0)."""
    alleles = parse_genotype(gt)
    if len(alleles) != 2 or alleles[0] is None or alleles[1] is None:
        return False
    return (alleles[0] == 0 and alleles[1] == 1) or (alleles[0] == 1 and alleles[1] == 0)


def is_hom_alt(gt: str) -> bool:
    """Check if genotype is homozygous alternate (1/1)."""
    alleles = parse_genotype(gt)
    if len(alleles) != 2 or alleles[0] is None or alleles[1] is None:
        return False
    return alleles[0] == 1 and alleles[1] == 1


def is_ref(gt: str) -> bool:
    """Check if genotype is homozygous reference (0/0) or haploid reference (0)."""
    alleles = parse_genotype(gt)
    if not alleles or any(a is None for a in alleles):
        return False
    return all(a == 0 for a in alleles)


def is_variant(gt: str) -> bool:
    """Check if genotype is a valid non-missing call containing an alternate allele."""
    if not gt or is_missing(gt) or is_ref(gt):
        return False
    alleles = parse_genotype(gt)
    if not alleles or any(a is None for a in alleles):
        return False
    return all(a in (0, 1) for a in alleles) and any(a == 1 for a in alleles)


def is_missing(gt: str) -> bool:
    """Check if genotype is missing or unknown."""
    if not gt:
        return True
    s = gt.strip()
    if s in (".", "./.", ".|.", "NA", "N/A", "nan", "NaN", "None", "null", "NULL", ""):
        return True
    if "." in s:
        return True
    alleles = parse_genotype(s)
    if not alleles:
        return True
    return any(a is None for a in alleles)


def get_allele_count(gt: str) -> int:
    """
    Get the count of alternate alleles in the genotype.

    Parameters
    ----------
    gt : str
        Genotype string

    Returns
    -------
    int
        Number of alternate alleles
    """
    if is_missing(gt) or not gt:
        return 0
    alleles = parse_genotype(gt)
    if not all(a in (0, 1) for a in alleles if a is not None):
        return 0
    return sum(1 for a in alleles if a == 1)


def could_be_de_novo(child_gt: str, father_gt: str, mother_gt: str) -> bool:
    """
    Check if a variant could be de novo in the child.

    Parameters
    ----------
    child_gt : str
        Child's genotype
    father_gt : str
        Father's genotype
    mother_gt : str
        Mother's genotype

    Returns
    -------
    bool
        True if the variant appears de novo in the child
    """
    # Child must have variant
    if not is_variant(child_gt):
        return False

    # Both parents must be reference
    return not (not is_ref(father_gt) or not is_ref(mother_gt))

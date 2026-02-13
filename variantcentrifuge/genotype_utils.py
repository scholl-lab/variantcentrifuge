"""
Genotype utility functions for inheritance analysis.

This module provides functions to parse and analyze genotype strings
commonly found in VCF files.
"""


def parse_genotype(gt: str) -> tuple[int | None, int | None]:
    """
    Parse a genotype string into allele indices.

    Parameters
    ----------
    gt : str
        Genotype string (e.g., "0/1", "1|1", "./.")

    Returns
    -------
    Tuple[Optional[int], Optional[int]]
        Tuple of (allele1, allele2) where None indicates missing
    """
    if not gt or gt == "./.":
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
    allele1, allele2 = parse_genotype(gt)
    if allele1 is None or allele2 is None:
        return False
    return (allele1 == 0 and allele2 == 1) or (allele1 == 1 and allele2 == 0)


def is_hom_alt(gt: str) -> bool:
    """Check if genotype is homozygous alternate (1/1)."""
    allele1, allele2 = parse_genotype(gt)
    if allele1 is None or allele2 is None:
        return False
    return allele1 == 1 and allele2 == 1


def is_ref(gt: str) -> bool:
    """Check if genotype is homozygous reference (0/0)."""
    allele1, allele2 = parse_genotype(gt)
    if allele1 is None or allele2 is None:
        return False
    return allele1 == 0 and allele2 == 0


def is_variant(gt: str) -> bool:
    """Check if genotype contains any alternate allele."""
    return is_het(gt) or is_hom_alt(gt)


def is_missing(gt: str) -> bool:
    """Check if genotype is missing or unknown."""
    if not gt:
        return True
    allele1, allele2 = parse_genotype(gt)
    return allele1 is None or allele2 is None


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
        Number of alternate alleles (0, 1, or 2)
    """
    if is_missing(gt) or is_ref(gt):
        return 0
    elif is_het(gt):
        return 1
    elif is_hom_alt(gt):
        return 2
    else:
        return 0


def is_phased(gt: str) -> bool:
    """Check if genotype is phased (uses | separator)."""
    return "|" in gt


def get_genotype_type(gt: str) -> str:
    """
    Get a string representation of the genotype type.

    Parameters
    ----------
    gt : str
        Genotype string

    Returns
    -------
    str
        One of: 'ref', 'het', 'hom_alt', 'missing'
    """
    if is_missing(gt):
        return "missing"
    elif is_ref(gt):
        return "ref"
    elif is_het(gt):
        return "het"
    elif is_hom_alt(gt):
        return "hom_alt"
    else:
        return "missing"


def is_mendelian_consistent(child_gt: str, father_gt: str, mother_gt: str) -> bool:
    """
    Check if a child's genotype is consistent with Mendelian inheritance.

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
        True if consistent with Mendelian inheritance
    """
    # Parse genotypes
    child_alleles = parse_genotype(child_gt)
    father_alleles = parse_genotype(father_gt)
    mother_alleles = parse_genotype(mother_gt)

    # If any genotype is missing, we can't determine consistency
    if None in child_alleles or None in father_alleles or None in mother_alleles:
        return True  # Assume consistent if we can't determine

    # Get possible alleles from each parent
    father_possible = set(father_alleles)
    mother_possible = set(mother_alleles)

    # Check if child's alleles could come from parents
    # One allele should come from father, one from mother
    for child_allele1 in [child_alleles[0], child_alleles[1]]:
        for child_allele2 in [child_alleles[0], child_alleles[1]]:
            if child_allele1 != child_allele2:
                continue
            if (child_alleles[0] in father_possible and child_alleles[1] in mother_possible) or (
                child_alleles[1] in father_possible and child_alleles[0] in mother_possible
            ):
                return True

    return False


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
    if not is_ref(father_gt) or not is_ref(mother_gt):
        return False

    return True


def merge_genotypes(genotypes: list[str]) -> str:
    """
    Merge multiple genotypes into a consensus genotype.

    Parameters
    ----------
    genotypes : List[str]
        List of genotype strings

    Returns
    -------
    str
        Merged genotype string
    """
    # Filter out missing genotypes
    valid_gts = [gt for gt in genotypes if not is_missing(gt)]

    if not valid_gts:
        return "./."

    # Count alleles
    allele_counts: dict[int, int] = {}
    for gt in valid_gts:
        allele1, allele2 = parse_genotype(gt)
        if allele1 is not None:
            allele_counts[allele1] = allele_counts.get(allele1, 0) + 1
        if allele2 is not None:
            allele_counts[allele2] = allele_counts.get(allele2, 0) + 1

    # Get the two most common alleles
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: x[1], reverse=True)

    if len(sorted_alleles) == 0:
        return "./."
    elif len(sorted_alleles) == 1:
        allele = sorted_alleles[0][0]
        return f"{allele}/{allele}"
    else:
        allele1, allele2 = sorted_alleles[0][0], sorted_alleles[1][0]
        return f"{min(allele1, allele2)}/{max(allele1, allele2)}"

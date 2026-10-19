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
        if "/" in gt or "|" in gt:
            sep = "|" if "|" in gt else "/"
            parts = [p.strip() for p in gt.split(sep)]
            if len(parts) == 2 and all(p.isalpha() for p in parts if p):
                return True
        return False
    return all(a in (0, 1) for a in alleles) and any(a == 1 for a in alleles)


def is_missing(gt: str) -> bool:
    """Check if genotype is missing or unknown."""
    if not gt:
        return True
    s = gt.strip()
    if s in (".", "./.", ".|.", "NA", "nan", "NaN", "None", "null", "NULL", ""):
        return True
    if "." in s:
        return True
    alleles = parse_genotype(s)
    if not alleles:
        return True
    if any(a is None for a in alleles):
        is_alpha_call = ("/" in s or "|" in s) and all(
            p.strip().isalpha() for p in s.replace("|", "/").split("/")
        )
        return not is_alpha_call
    return False


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
        One of: 'ref', 'het', 'hom_alt', 'haploid_alt', 'missing'
    """
    if is_missing(gt):
        return "missing"
    elif is_ref(gt):
        return "ref"
    elif is_het(gt):
        return "het"
    elif is_hom_alt(gt):
        return "hom_alt"
    elif is_variant(gt):
        return "haploid_alt"
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
    if len(child_alleles) == 1:
        return child_alleles[0] in father_possible or child_alleles[0] in mother_possible

    if len(child_alleles) == 2:
        return (child_alleles[0] in father_possible and child_alleles[1] in mother_possible) or (
            child_alleles[1] in father_possible and child_alleles[0] in mother_possible
        )

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
    return not (not is_ref(father_gt) or not is_ref(mother_gt))


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
        if any(gt.strip() == "." for gt in genotypes if isinstance(gt, str)) and not any(
            gt.strip() in ("./.", ".|.") for gt in genotypes if isinstance(gt, str)
        ):
            return "."
        return "./."

    valid_alleles = [parse_genotype(gt) for gt in valid_gts]
    is_all_haploid = all(len(alleles) == 1 for alleles in valid_alleles)

    # Count alleles
    allele_counts: dict[int, int] = {}
    for alleles in valid_alleles:
        for allele in alleles:
            if allele is not None:
                allele_counts[allele] = allele_counts.get(allele, 0) + 1

    # Get the most common alleles
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: x[1], reverse=True)

    if len(sorted_alleles) == 0:
        return "." if is_all_haploid else "./."

    if is_all_haploid:
        if len(sorted_alleles) == 1 or sorted_alleles[0][1] > sorted_alleles[1][1]:
            return str(sorted_alleles[0][0])
        return "."

    if len(sorted_alleles) == 1:
        allele = sorted_alleles[0][0]
        return f"{allele}/{allele}"
    allele1, allele2 = sorted_alleles[0][0], sorted_alleles[1][0]
    return f"{min(allele1, allele2)}/{max(allele1, allele2)}"

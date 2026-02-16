"""
Synthetic data generators for performance benchmarking.

All generated data is fully synthetic with reproducible seeding. No real patient
identifiers or cohort names are used. All samples are named SAMPLE_NNNN, genes
are GENE_NNNN, and pedigree members use CHILD/FATHER/MOTHER prefixes.
"""

import numpy as np
import pandas as pd


def generate_synthetic_variants(n_variants: int, n_samples: int, seed: int = 42) -> pd.DataFrame:
    """
    Generate synthetic variant DataFrame matching variantcentrifuge structure.

    Creates realistic genomic variant data with proper allele frequency distributions,
    effect types, and impact levels. All identifiers are fully synthetic.

    Parameters
    ----------
    n_variants : int
        Number of variants to generate
    n_samples : int
        Number of samples to generate genotypes for
    seed : int, default=42
        Random seed for reproducibility

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: CHROM, POS, REF, ALT, GT, GENE, FILTER, EFFECT, IMPACT
        Sorted by CHROM (natural order) then POS

    Examples
    --------
    >>> df = generate_synthetic_variants(100, 10, seed=42)
    >>> df.shape
    (100, 9)
    >>> list(df.columns)
    ['CHROM', 'POS', 'REF', 'ALT', 'GT', 'GENE', 'FILTER', 'EFFECT', 'IMPACT']
    >>> df1 = generate_synthetic_variants(50, 5, seed=123)
    >>> df2 = generate_synthetic_variants(50, 5, seed=123)
    >>> df1.equals(df2)
    True
    """
    rng = np.random.default_rng(seed)

    # Generate chromosomes with realistic distribution (favor chr1-10)
    chrom_choices = [str(i) for i in range(1, 23)] + ["X", "Y"]
    # Create weights that sum to exactly 1.0: chr1 gets extra to balance
    chrom_weights = np.array([0.135] + [0.075] * 9 + [0.015] * 12 + [0.005, 0.005])
    # Normalize to ensure exact sum of 1.0 (handles floating point precision)
    chrom_weights = chrom_weights / chrom_weights.sum()
    chroms = rng.choice(chrom_choices, size=n_variants, p=chrom_weights)

    # Generate positions (1 to 250M)
    positions = rng.integers(1, 250_000_000, size=n_variants)

    # Generate REF/ALT nucleotides (ensure REF != ALT)
    bases = ["A", "C", "G", "T"]
    refs = rng.choice(bases, size=n_variants)
    alts = np.array([rng.choice([b for b in bases if b != ref]) for ref in refs])

    # Generate genes (~n_variants/10 unique genes so avg ~10 variants per gene)
    n_genes = max(1, n_variants // 10)
    gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
    genes = rng.choice(gene_names, size=n_variants)

    # Generate FILTER (95% PASS, 5% LowQual)
    filters = rng.choice(["PASS", "LowQual"], size=n_variants, p=[0.95, 0.05])

    # Generate EFFECT with realistic weights
    effect_choices = [
        "missense_variant",
        "synonymous_variant",
        "frameshift_variant",
        "stop_gained",
        "splice_region_variant",
    ]
    effect_weights = [0.40, 0.30, 0.10, 0.05, 0.15]
    effects = rng.choice(effect_choices, size=n_variants, p=effect_weights)

    # Map EFFECT -> IMPACT deterministically
    effect_to_impact = {
        "frameshift_variant": "HIGH",
        "stop_gained": "HIGH",
        "missense_variant": "MODERATE",
        "synonymous_variant": "LOW",
        "splice_region_variant": "MODIFIER",
    }
    impacts = np.array([effect_to_impact[e] for e in effects])

    # Generate genotypes with realistic allele frequency distribution
    # 90% rare (95/4/1 0-0/0-1/1-1), 8% uncommon (80/15/5), 2% common (50/40/10)
    genotype_choices = ["0/0", "0/1", "1/1"]
    genotypes_list = []

    for _ in range(n_variants):
        variant_class = rng.choice(["rare", "uncommon", "common"], p=[0.90, 0.08, 0.02])

        if variant_class == "rare":
            gt_weights = [0.95, 0.04, 0.01]
        elif variant_class == "uncommon":
            gt_weights = [0.80, 0.15, 0.05]
        else:  # common
            gt_weights = [0.50, 0.40, 0.10]

        sample_gts = rng.choice(genotype_choices, size=n_samples, p=gt_weights)
        genotypes_list.append(",".join(sample_gts))

    # Create DataFrame
    df = pd.DataFrame(
        {
            "CHROM": chroms,
            "POS": positions,
            "REF": refs,
            "ALT": alts,
            "GT": genotypes_list,
            "GENE": genes,
            "FILTER": filters,
            "EFFECT": effects,
            "IMPACT": impacts,
        }
    )

    # Sort by CHROM (natural order) then POS
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order["X"] = 23
    chrom_order["Y"] = 24
    df["_chrom_sort"] = df["CHROM"].map(chrom_order)
    df = df.sort_values(["_chrom_sort", "POS"]).drop(columns=["_chrom_sort"]).reset_index(drop=True)

    return df


def generate_synthetic_pedigree(n_samples: int = 3, seed: int = 42) -> dict:
    """
    Generate synthetic pedigree matching variantcentrifuge format.

    Creates a trio (child, father, mother) plus additional unrelated samples if requested.
    All identifiers are fully synthetic.

    Parameters
    ----------
    n_samples : int, default=3
        Total number of samples (minimum 3 for trio)
    seed : int, default=42
        Random seed for reproducibility

    Returns
    -------
    dict
        Pedigree data in format {sample_id: {sample_id, father_id, mother_id, sex, affected_status}}

    Examples
    --------
    >>> ped = generate_synthetic_pedigree(5, seed=42)
    >>> len(ped)
    5
    >>> 'CHILD_001' in ped
    True
    >>> ped['CHILD_001']['affected_status']
    '2'
    """
    if n_samples < 3:
        raise ValueError("Minimum 3 samples required for a trio")

    rng = np.random.default_rng(seed)

    pedigree = {}

    # Core trio
    pedigree["CHILD_001"] = {
        "sample_id": "CHILD_001",
        "father_id": "FATHER_001",
        "mother_id": "MOTHER_001",
        "sex": "1" if rng.random() < 0.5 else "2",  # 1=male, 2=female
        "affected_status": "2",  # 2=affected
    }

    pedigree["FATHER_001"] = {
        "sample_id": "FATHER_001",
        "father_id": "0",
        "mother_id": "0",
        "sex": "1",  # male
        "affected_status": "1",  # 1=unaffected
    }

    pedigree["MOTHER_001"] = {
        "sample_id": "MOTHER_001",
        "father_id": "0",
        "mother_id": "0",
        "sex": "2",  # female
        "affected_status": "1",  # 1=unaffected
    }

    # Additional samples as unrelated individuals
    for i in range(4, n_samples + 1):
        sample_id = f"SAMPLE_{i - 3:04d}"  # SAMPLE_0001, SAMPLE_0002, etc.
        pedigree[sample_id] = {
            "sample_id": sample_id,
            "father_id": "0",
            "mother_id": "0",
            "sex": "1" if rng.random() < 0.5 else "2",
            "affected_status": "1" if rng.random() < 0.8 else "2",  # 80% unaffected
        }

    return pedigree


def generate_synthetic_scoring_config() -> dict:
    """
    Generate synthetic scoring configuration for benchmark tests.

    Returns minimal scoring config matching variantcentrifuge format with
    variable mappings and formulas compatible with pd.eval.

    Returns
    -------
    dict
        Scoring configuration with 'variables' and 'formulas' keys

    Examples
    --------
    >>> config = generate_synthetic_scoring_config()
    >>> 'variables' in config and 'formulas' in config
    True
    """
    return {
        "variables": {
            "impact_var": {
                "IMPACT": {
                    "HIGH": 3,
                    "MODERATE": 2,
                    "LOW": 1,
                    "MODIFIER": 0,
                }
            },
            "effect_var": {
                "EFFECT": {
                    "stop_gained": 5,
                    "frameshift_variant": 4,
                    "missense_variant": 2,
                    "splice_region_variant": 1,
                    "synonymous_variant": 0,
                }
            },
        },
        "formulas": [
            {"impact_score": "impact_var"},
            {"combined_score": "impact_var * 2 + effect_var"},
        ],
    }


def generate_synthetic_gene_burden_data(
    n_variants: int, n_samples: int, n_genes: int = 50, seed: int = 42
) -> tuple[pd.DataFrame, set, set]:
    """
    Generate synthetic data for gene burden analysis benchmarks.

    Creates DataFrame with GT column in gene_burden format:
    "SAMPLE_0001(0/1);SAMPLE_0002(0/0);..."

    Parameters
    ----------
    n_variants : int
        Number of variants to generate
    n_samples : int
        Number of samples (split 50/50 into cases and controls)
    n_genes : int, default=50
        Number of unique genes to distribute variants across
    seed : int, default=42
        Random seed for reproducibility

    Returns
    -------
    tuple
        (df, case_samples, control_samples) where df has GENE and GT columns,
        and sample sets are disjoint

    Examples
    --------
    >>> df, cases, controls = generate_synthetic_gene_burden_data(100, 20, n_genes=10, seed=42)
    >>> df.shape[0]
    100
    >>> len(cases) + len(controls)
    20
    >>> len(cases & controls)
    0
    """
    rng = np.random.default_rng(seed)

    # Generate gene names
    gene_names = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]
    genes = rng.choice(gene_names, size=n_variants)

    # Split samples into cases and controls (50/50)
    n_cases = n_samples // 2
    n_controls = n_samples - n_cases

    case_samples = {f"SAMPLE_{i:04d}" for i in range(1, n_cases + 1)}
    control_samples = {f"SAMPLE_{i:04d}" for i in range(n_cases + 1, n_samples + 1)}

    # Generate genotypes in gene_burden format
    genotype_choices = ["0/0", "0/1", "1/1"]
    genotypes_list = []

    for _ in range(n_variants):
        # Case samples: higher variant frequency
        case_gts = rng.choice(genotype_choices, size=n_cases, p=[0.70, 0.25, 0.05])
        # Control samples: lower variant frequency
        control_gts = rng.choice(genotype_choices, size=n_controls, p=[0.90, 0.08, 0.02])

        # Combine into gene_burden GT format
        all_samples = list(case_samples) + list(control_samples)
        all_gts = list(case_gts) + list(control_gts)
        gt_str = ";".join(f"{s}({gt})" for s, gt in zip(all_samples, all_gts, strict=False))
        genotypes_list.append(gt_str)

    # Create DataFrame with required columns for gene burden analysis
    # The perform_gene_burden_analysis function expects these columns:
    # - GENE: gene identifier
    # - GT: genotypes in format "SAMPLE_0001(0/1);SAMPLE_0002(0/0);..."
    # - proband_count: total number of case samples
    # - control_count: total number of control samples
    # Note: The function will parse GT and aggregate per gene internally
    df = pd.DataFrame(
        {
            "GENE": genes,
            "GT": genotypes_list,
            "proband_count": n_cases,
            "control_count": n_controls,
        }
    )

    return df, case_samples, control_samples

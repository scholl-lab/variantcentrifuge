# File: variantcentrifuge/helpers.py
# Location: variantcentrifuge/variantcentrifuge/helpers.py

"""
Helper functions for analyze_variants and related processes.

Provides:
- Case/control assignment
- Phenotype map building
- Genotype parsing and allele count conversion
- Sample and phenotype classification logic
"""

import logging
import sys
from collections import defaultdict
from typing import Dict, Any, Set, Tuple, List
import pandas as pd

logger = logging.getLogger("variantcentrifuge")


def determine_case_control_sets(all_samples: Set[str], cfg: Dict[str, Any], df: pd.DataFrame) -> Tuple[Set[str], Set[str]]:
    """
    Determine case/control sample sets based on configuration.

    Logic:
    - If explicit case/control samples are provided, use them directly.
    - Else if phenotype terms are given, classify samples based on those.
    - Else, all samples become controls.

    Parameters
    ----------
    all_samples : set of str
        The full set of sample names.
    cfg : dict
        Configuration dictionary that may include "case_samples", "control_samples",
        "case_phenotypes", "control_phenotypes".
    df : pd.DataFrame
        DataFrame with variant information, potentially containing phenotypes.

    Returns
    -------
    Tuple[set, set]
        (set_of_case_samples, set_of_control_samples)

    Raises
    ------
    SystemExit
        If no valid configuration for assigning case/control sets can be determined.
    """
    logger.debug("Determining case/control sets...")

    case_samples = set(cfg.get("case_samples", []))
    control_samples = set(cfg.get("control_samples", []))

    # Step 1: If explicit sets are provided
    if case_samples or control_samples:
        logger.debug("Explicit sample sets provided: %d cases, %d controls", len(case_samples), len(control_samples))
        if case_samples and not control_samples:
            control_samples = all_samples - case_samples
        elif control_samples and not case_samples:
            case_samples = all_samples - control_samples
        return case_samples, control_samples

    # Step 2: Phenotype-based logic
    case_terms = cfg.get("case_phenotypes", [])
    control_terms = cfg.get("control_phenotypes", [])

    if not case_terms and not control_terms:
        # Step 3: No criteria, all controls
        logger.debug("No phenotype terms or sets provided. All samples = controls.")
        return set(), all_samples

    sample_phenotype_map = build_sample_phenotype_map(df)
    logger.debug("Phenotype map has %d samples.", len(sample_phenotype_map))

    classified_cases = set()
    classified_controls = set()

    for s in all_samples:
        phenos = sample_phenotype_map.get(s, set())
        match_case = any(p in phenos for p in case_terms) if case_terms else False
        match_control = any(p in phenos for p in control_terms) if control_terms else False

        if case_terms and not control_terms:
            # Only case terms
            if match_case:
                classified_cases.add(s)
            else:
                classified_controls.add(s)
        elif control_terms and not case_terms:
            # Only control terms
            if match_control:
                classified_controls.add(s)
            else:
                classified_cases.add(s)
        else:
            # Both sets of terms
            if match_case and not match_control:
                classified_cases.add(s)
            elif match_control and not match_case:
                classified_controls.add(s)
            # If matches both or none, sample is not classified

    if case_terms and len(classified_cases) == 0:
        logger.warning("No samples match the case phenotype terms.")
    if control_terms and len(classified_controls) == 0:
        logger.warning("No samples match the control phenotype terms.")

    logger.debug("Classified %d cases, %d controls.", len(classified_cases), len(classified_controls))
    return classified_cases, classified_controls


def build_sample_phenotype_map(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """
    Build a map of sample -> set_of_phenotypes from the 'phenotypes' column if present.

    If multiple phenotype groups appear and multiple samples in GT are present:
    - If counts match, assign phenotypes group-wise.
    - Otherwise, assign all phenotypes to all samples in that variant line.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing a 'phenotypes' column and 'GT' column.

    Returns
    -------
    dict of {str: set of str}
        A dictionary mapping sample names to a set of phenotypes.
    """
    if "phenotypes" not in df.columns:
        return {}

    from collections import defaultdict
    sample_phenos = defaultdict(set)

    for idx, row in df.iterrows():
        pheno_str = row.get("phenotypes", "")
        if not isinstance(pheno_str, str) or not pheno_str.strip():
            continue

        gt_val = row.get("GT", "")
        sample_entries = gt_val.split(";") if gt_val and isinstance(gt_val, str) else []
        sample_names = []
        for entry in sample_entries:
            entry = entry.strip()
            if not entry:
                continue
            sname, _ = extract_sample_and_genotype(entry)
            if sname:
                sample_names.append(sname)

        pheno_groups = [pg.strip() for pg in pheno_str.split(";") if pg.strip()]

        if len(pheno_groups) > 1:
            # Multiple phenotype groups
            if len(pheno_groups) == len(sample_names):
                # Perfect alignment
                for sname, pgroup in zip(sample_names, pheno_groups):
                    phenos = {p.strip() for p in pgroup.split(",") if p.strip()}
                    sample_phenos[sname].update(phenos)
            else:
                # Mismatch in counts, assign all to all samples
                logger.warning(
                    "Phenotype groups (%d) != sample count (%d) at row %d. "
                    "Assigning all phenotypes to all samples in this row.",
                    len(pheno_groups),
                    len(sample_names),
                    idx
                )
                combined_phenos = set()
                for pgroup in pheno_groups:
                    combined_phenos.update({p.strip() for p in pgroup.split(",") if p.strip()})
                for sname in sample_names:
                    sample_phenos[sname].update(combined_phenos)
        else:
            # Only one phenotype group (or none)
            phenos = {p.strip() for p in pheno_str.split(",") if p.strip()}
            for sname in sample_names:
                sample_phenos[sname].update(phenos)

    return dict(sample_phenos)


def assign_case_control_counts(df: pd.DataFrame, case_samples: Set[str], control_samples: Set[str],
                               all_samples: Set[str]) -> pd.DataFrame:
    """
    Assign case/control counts, allele counts, and homozygous variant counts per variant.

    Creates columns:
    - proband_count/control_count: total number of case/control samples
    - proband_variant_count/control_variant_count: number of case/control samples with a variant allele
    - proband_allele_count/control_allele_count: sum of variant alleles in case/control samples
    - proband_homozygous_count/control_homozygous_count: how many case/control samples have a homozygous variant (1/1)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of variants with a "GT" column listing variants per sample.
    case_samples : set of str
        Set of samples classified as cases.
    control_samples : set of str
        Set of samples classified as controls.
    all_samples : set of str
        All samples present in the VCF.

    Returns
    -------
    pd.DataFrame
        DataFrame with assigned case/control counts and alleles, including homozygous counts.
    """
    logger.debug("Assigning case/control counts to variants...")

    total_proband = len(case_samples)
    total_control = len(control_samples)
    logger.debug("Total proband samples: %d, Total control samples: %d", total_proband, total_control)

    proband_variant_count_list = []
    control_variant_count_list = []
    proband_allele_count_list = []
    control_allele_count_list = []
    p_hom_count_list = []
    c_hom_count_list = []

    for idx, val in enumerate(df["GT"]):
        # Log progress every 1000 variants
        if idx % 1000 == 0:
            logger.debug("Processing variant row %d for allele counts...", idx)

        samples_with_variant = {}
        if isinstance(val, str) and val.strip():
            for s in val.split(";"):
                s = s.strip()
                if not s:
                    continue
                sample_name, genotype = extract_sample_and_genotype(s)
                if sample_name is not None:
                    samples_with_variant[sample_name] = genotype

        p_variant_count = 0
        c_variant_count = 0
        p_allele_count = 0
        c_allele_count = 0
        p_hom_count = 0
        c_hom_count = 0

        # For each sample, if not in samples_with_variant, assume genotype=0/0
        for sample_name in all_samples:
            genotype = samples_with_variant.get(sample_name, "0/0")
            allele_count = genotype_to_allele_count(genotype)

            # Check if homozygous variant
            is_hom_variant = (genotype == "1/1")

            if sample_name in case_samples:
                if allele_count > 0:
                    p_variant_count += 1
                    p_allele_count += allele_count
                    if is_hom_variant:
                        p_hom_count += 1
            elif sample_name in control_samples:
                if allele_count > 0:
                    c_variant_count += 1
                    c_allele_count += allele_count
                    if is_hom_variant:
                        c_hom_count += 1

        proband_variant_count_list.append(p_variant_count)
        control_variant_count_list.append(c_variant_count)
        proband_allele_count_list.append(p_allele_count)
        control_allele_count_list.append(c_allele_count)
        p_hom_count_list.append(p_hom_count)
        c_hom_count_list.append(c_hom_count)

    df["proband_count"] = total_proband
    df["control_count"] = total_control
    df["proband_variant_count"] = proband_variant_count_list
    df["control_variant_count"] = control_variant_count_list
    df["proband_allele_count"] = proband_allele_count_list
    df["control_allele_count"] = control_allele_count_list
    df["proband_homozygous_count"] = p_hom_count_list
    df["control_homozygous_count"] = c_hom_count_list

    logger.debug("Case/control counts assigned.")
    if len(df) > 0:
        logger.debug("Example for first row: %s", df.iloc[0][[
            "proband_count", "proband_variant_count", "proband_allele_count",
            "proband_homozygous_count",
            "control_count", "control_variant_count", "control_allele_count",
            "control_homozygous_count"
        ]].to_dict())

    return df


def extract_sample_and_genotype(sample_field: str) -> Tuple[str, str]:
    """
    Extract sample name and genotype from a field like 'sample(0/1)'.

    If parentheses are missing, assume no genotype is specified -> no variant (0/0).

    Parameters
    ----------
    sample_field : str
        A string like 'sample(0/1)' or 'sample'.

    Returns
    -------
    (str, str)
        (sample_name, genotype)
    """
    start_idx = sample_field.find("(")
    end_idx = sample_field.find(")")
    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        sample_name = sample_field[:start_idx].strip()
        genotype = sample_field[start_idx + 1:end_idx].strip()
        return sample_name, genotype
    else:
        return sample_field.strip(), ""


def genotype_to_allele_count(genotype: str) -> int:
    """
    Convert genotype string to allele count:
    - '1/1' -> 2
    - '0/1' or '1/0' -> 1
    - '0/0' or '' -> 0

    Parameters
    ----------
    genotype : str
        Genotype string, expected to be one of '0/0', '0/1', '1/0', '1/1', or ''.

    Returns
    -------
    int
        The allele count for the given genotype.
    """
    if genotype == "1/1":
        return 2
    elif genotype in ["0/1", "1/0"]:
        return 1
    return 0

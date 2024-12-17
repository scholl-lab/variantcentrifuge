"""
Helper functions for analyze_variants and related processes.

Provides:
- Case/control assignment
- Phenotype map building
- Genotype parsing and allele count conversion
- Sample and phenotype classification logic
"""

import logging
from collections import defaultdict

logger = logging.getLogger("variantcentrifuge")


def determine_case_control_sets(all_samples, cfg, df):
    """
    Determine case/control sample sets based on configuration.

    Logic:
    - If explicit case/control samples provided, use them directly.
    - Else if phenotype terms given, classify based on those.
    - Else, all samples become controls.

    Returns
    -------
    (set_of_case_samples, set_of_control_samples)
    """
    logger.debug("Determining case/control sets...")

    case_samples = set(cfg.get("case_samples", []))
    control_samples = set(cfg.get("control_samples", []))

    # Step 1: Explicit sets
    if case_samples or control_samples:
        logger.debug(f"Explicit sample sets: {len(case_samples)} cases, {len(control_samples)} controls")
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
        logger.debug("No phenotype or sets given. All samples = controls.")
        return set(), all_samples

    sample_phenotype_map = build_sample_phenotype_map(df)
    logger.debug(f"Phenotype map with {len(sample_phenotype_map)} samples.")

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
            # Else excluded

    if case_terms and len(classified_cases) == 0:
        logger.warning("No samples match the case phenotype terms.")
    if control_terms and len(classified_controls) == 0:
        logger.warning("No samples match the control phenotype terms.")

    logger.debug(f"Classified {len(classified_cases)} cases, {len(classified_controls)} controls.")
    return classified_cases, classified_controls


def build_sample_phenotype_map(df):
    """
    Build a map of sample -> set_of_phenotypes from the 'phenotypes' column if present.

    If multiple phenotype groups and multiple samples in GT, align them if counts match.
    Otherwise, assign all phenotypes to all samples in that variant line.

    Returns
    -------
    dict
        {sample_name: set_of_phenotypes}
    """
    if "phenotypes" not in df.columns:
        return {}

    sample_phenos = defaultdict(set)

    for idx, row in df.iterrows():
        pheno_str = row.get("phenotypes", "")
        if not isinstance(pheno_str, str) or not pheno_str.strip():
            continue

        gt_val = row["GT"]
        sample_entries = gt_val.split(";") if gt_val and isinstance(gt_val, str) else []
        sample_names = []
        for entry in sample_entries:
            entry = entry.strip()
            if not entry:
                continue
            sname = entry
            paren_idx = sname.find("(")
            if paren_idx != -1:
                sname = sname[:paren_idx].strip()
            if sname:
                sample_names.append(sname)

        pheno_groups = pheno_str.split(";")
        pheno_groups = [pg.strip() for pg in pheno_groups if pg.strip()]

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
                    f"Phenotype groups ({len(pheno_groups)}) != sample count ({len(sample_names)}) at row {idx}. "
                    "Assigning all phenotypes to all samples in this row."
                )
                combined_phenos = set()
                for pgroup in pheno_groups:
                    combined_phenos.update({p.strip() for p in pgroup.split(",") if p.strip()})
                for sname in sample_names:
                    sample_phenos[sname].update(combined_phenos)
        else:
            # Only one phenotype group or none
            phenos = {p.strip() for p in pheno_str.split(",") if p.strip()}
            for sname in sample_names:
                sample_phenos[sname].update(phenos)

    return sample_phenos


def assign_case_control_counts(df, case_samples, control_samples, all_samples):
    """
    Assign case/control counts and allele counts per variant.

    Creates columns:
    - proband_count/control_count: total number of case/control samples
    - proband_variant_count/control_variant_count: how many case/control samples have the variant
    - proband_allele_count/control_allele_count: sum of variant alleles in case/control samples

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of variants with a GT column listing variants per sample.
    case_samples : set
        Set of samples classified as cases.
    control_samples : set
        Set of samples classified as controls.
    all_samples : set
        All samples present in the VCF.

    Returns
    -------
    pd.DataFrame
        DataFrame with assigned case/control counts and alleles.
    """
    logger.debug("Assigning case/control counts to variants...")

    total_proband = len(case_samples)
    total_control = len(control_samples)
    logger.debug(f"Total proband samples: {total_proband}, Total control samples: {total_control}")

    proband_variant_count_list = []
    control_variant_count_list = []
    proband_allele_count_list = []
    control_allele_count_list = []

    for idx, val in enumerate(df["GT"]):
        # Log progress every 1000 variants
        if idx % 1000 == 0:
            logger.debug(f"Processing variant row {idx} for allele counts...")

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

        # For each sample, if not in samples_with_variant, assume genotype=0/0
        for sample_name in all_samples:
            genotype = samples_with_variant.get(sample_name, "0/0")
            allele_count = genotype_to_allele_count(genotype)

            if sample_name in case_samples:
                if allele_count > 0:
                    p_variant_count += 1
                    p_allele_count += allele_count
            elif sample_name in control_samples:
                if allele_count > 0:
                    c_variant_count += 1
                    c_allele_count += allele_count

        proband_variant_count_list.append(p_variant_count)
        control_variant_count_list.append(c_variant_count)
        proband_allele_count_list.append(p_allele_count)
        control_allele_count_list.append(c_allele_count)

    df["proband_count"] = total_proband
    df["control_count"] = total_control
    df["proband_variant_count"] = proband_variant_count_list
    df["control_variant_count"] = control_variant_count_list
    df["proband_allele_count"] = proband_allele_count_list
    df["control_allele_count"] = control_allele_count_list

    logger.debug("Case/control counts assigned. Example for first row:")
    if len(df) > 0:
        logger.debug(df.iloc[0][[
            "proband_count", "proband_variant_count", "proband_allele_count",
            "control_count", "control_variant_count", "control_allele_count"
        ]].to_dict())

    return df


def extract_sample_and_genotype(sample_field):
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
        genotype = sample_field[start_idx+1:end_idx].strip()
        return sample_name, genotype
    else:
        return sample_field.strip(), ""


def genotype_to_allele_count(genotype):
    """
    Convert genotype string to allele count:
    - '1/1' -> 2
    - '0/1' or '1/0' -> 1
    - '0/0' or '' -> 0
    """
    if genotype == "1/1":
        return 2
    elif genotype in ["0/1", "1/0"]:
        return 1
    return 0

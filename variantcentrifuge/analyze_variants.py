# File: variantcentrifuge/analyze_variants.py
# Location: variantcentrifuge/variantcentrifuge/analyze_variants.py

"""
Variant analysis module for gene burden and other statistics.

This module:
- Reads a TSV of variants and their annotations (including a GT column listing sample genotypes).
- Classifies samples into case/control sets based on user input (sample lists or phenotype terms).
- If no case/control criteria are provided, defaults to making all samples controls.
- Computes variant statistics and optionally performs a gene burden analysis via Fisher's exact test.

Fallback Behavior:
- If case or control groups are not provided by any method (phenotypes, sample lists), 
  all samples are considered controls by default.

Added Debug Logging:
- Logs at each major logical step and data processing point.
- Logs the size of sample sets and the outcome of classification.

Added Comments:
- Each function now has a docstring and inline comments for clarity.
"""

import io
import pandas as pd
from collections import defaultdict
import logging
from math import isnan

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

logger = logging.getLogger("variantcentrifuge")

def analyze_variants(lines, cfg):
    """
    Analyze variants, perform gene burden analysis if requested, and add debug logging.

    Steps:
    - Parse TSV input into a DataFrame.
    - Retrieve the full sample list from cfg["sample_list"] (extracted from the VCF header).
    - Determine case/control sets based on cfg:
      * If case_samples/control_samples given, use them directly.
      * Else if case_phenotypes/control_phenotypes given, classify samples accordingly.
      * Else, all samples become controls.
    - Assign allele counts per variant (row) into proband_*/control_* columns.
    - Compute and write basic stats, optionally comprehensive stats.
    - If requested, perform Fisher's exact test for gene burden.

    Parameters
    ----------
    lines : iterator of str
        Input lines representing a TSV with variant data.
    cfg : dict
        Configuration dictionary with keys:
        - sample_list: comma-separated full sample list from VCF
        - case_samples, control_samples (list of strings)
        - case_phenotypes, control_phenotypes (list of strings)
        - perform_gene_burden (bool)
        - no_stats (bool)
        - stats_output_file (str)

    Yields
    ------
    str
        Processed lines of output TSV or gene burden results.
    """
    perform_gene_burden = cfg.get("perform_gene_burden", False)
    stats_output_file = cfg.get("stats_output_file")
    no_stats = cfg.get("no_stats", False)

    logger.debug(f"analyze_variants: perform_gene_burden={perform_gene_burden}, "
                 f"stats_output_file={stats_output_file}, no_stats={no_stats}")

    # Read all input lines
    text_data = "".join(line for line in lines)
    if not text_data.strip():
        logger.debug("No input data provided to analyze_variants.")
        return

    df = pd.read_csv(io.StringIO(text_data), sep="\t", dtype=str)
    if len(df) == 0:
        logger.debug("Empty DataFrame after reading input in analyze_variants.")
        return

    required_columns = ["CHROM", "POS", "REF", "ALT", "GENE", "GT"]
    missing_columns = [c for c in required_columns if c not in df.columns]
    if missing_columns:
        logger.error(f"Missing required columns: {', '.join(missing_columns)}. Returning unchanged lines.")
        for line in text_data.strip().split("\n"):
            yield line
        return

    # Ensure we have the full sample list from the VCF
    if "sample_list" not in cfg or not cfg["sample_list"].strip():
        logger.error("No sample_list found in cfg. Unable to determine the full sample set.")
        for line in text_data.strip().split("\n"):
            yield line
        return

    # all_samples is the full set of samples from the VCF header
    all_samples = set(cfg["sample_list"].split(","))
    logger.debug(f"Total samples from VCF header: {len(all_samples)}")

    # Determine case/control sets
    case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)
    logger.debug(f"Number of case samples: {len(case_samples)}; Number of control samples: {len(control_samples)}")

    # Assign case/control counts to the variants
    df = assign_case_control_counts(df, case_samples, control_samples, all_samples)

    # Compute basic stats
    logger.debug("Computing basic statistics...")
    basic_stats_df = compute_basic_stats(df, all_samples)
    logger.debug("Basic statistics computed.")

    # Write basic stats if requested
    if stats_output_file:
        logger.debug(f"Writing basic stats to {stats_output_file}")
        basic_stats_df.to_csv(stats_output_file, sep="\t", index=False, header=True, mode="w")

    # Compute comprehensive stats if requested
    if not no_stats:
        logger.debug("Computing comprehensive gene-level statistics...")
        comprehensive_df = compute_gene_stats(df)
        impact_summary = compute_impact_summary(df)
        variant_type_summary = compute_variant_type_summary(df)
        combined_stats = merge_and_format_stats(comprehensive_df, impact_summary, variant_type_summary)
        logger.debug("Comprehensive statistics computed.")

        comp_stats_list = []
        for _, row in combined_stats.iterrows():
            gene = row["GENE"]
            for col in combined_stats.columns:
                if col == "GENE":
                    continue
                val = row[col]
                comp_stats_list.append([f"{gene}_{col}", str(val)])
        if comp_stats_list and stats_output_file:
            logger.debug("Appending comprehensive stats to the stats file.")
            comp_stats_df = pd.DataFrame(comp_stats_list, columns=["metric", "value"])
            comp_stats_df.to_csv(stats_output_file, sep="\t", index=False, header=False, mode="a")
    else:
        logger.debug("No comprehensive stats requested (no_stats=True).")

    # Perform gene burden if requested
    if perform_gene_burden:
        logger.debug("Performing gene burden analysis (Fisher's exact test)...")
        if fisher_exact is None:
            logger.error("scipy not available for Fisher test, cannot perform gene burden.")
            for line in text_data.strip().split("\n"):
                yield line
            return
        grouped = df.groupby("GENE").apply(gene_burden_fisher)
        burden_text = grouped.to_csv(sep="\t", index=False)
        logger.debug("Gene burden analysis complete.")
        for line in burden_text.strip().split("\n"):
            yield line
    else:
        logger.debug("Gene burden not requested, returning processed data.")
        out_str = df.to_csv(sep="\t", index=False)
        for line in out_str.strip().split("\n"):
            yield line

def determine_case_control_sets(all_samples, cfg, df):
    """
    Determine which samples are cases and which are controls.

    Priority:
    1. If case_samples or control_samples explicitly provided, use them directly.
       - If only case_samples are provided: non-case are controls.
       - If only control_samples are provided: non-control are cases.
       - If both provided: exclude samples not in either set.
    2. Else, use phenotype-based logic if phenotype terms are provided.
       - If only case terms: samples matching any case term -> case; others -> control
       - If only control terms: samples matching any control term -> control; others -> case
       - If both sets of terms: 
         * match case terms -> case, 
         * match control terms -> control, 
         * match neither -> excluded
    3. If no phenotype or sample sets are provided, all samples become controls by default.

    Returns
    -------
    (set_of_case_samples, set_of_control_samples)
    """

    logger.debug("Determining case/control sets...")

    case_samples = set(cfg.get("case_samples", []))
    control_samples = set(cfg.get("control_samples", []))

    # Step 1: Explicit sample sets
    if case_samples or control_samples:
        logger.debug(f"Explicit sample sets provided: {len(case_samples)} cases, {len(control_samples)} controls")
        if case_samples and not control_samples:
            logger.debug("Only case samples given, assigning all others as controls.")
            control_samples = all_samples - case_samples
        elif control_samples and not case_samples:
            logger.debug("Only control samples given, assigning all others as cases.")
            case_samples = all_samples - control_samples
        # If both sets provided, just use them as is (excluded are those not in either)
        return case_samples, control_samples

    # Step 2: Phenotype-based logic
    case_terms = cfg.get("case_phenotypes", [])
    control_terms = cfg.get("control_phenotypes", [])

    if not case_terms and not control_terms:
        # Step 3: No sets or phenotypes given, all samples become controls
        logger.debug("No phenotype or explicit sets given. All samples are controls.")
        return set(), all_samples

    # Build phenotype map
    sample_phenotype_map = build_sample_phenotype_map(df)
    logger.debug(f"Phenotype map built with {len(sample_phenotype_map)} samples having phenotypes.")

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
            # Both case and control terms
            if match_case and not match_control:
                classified_cases.add(s)
            elif match_control and not match_case:
                classified_controls.add(s)
            # Neither matches -> excluded

    if case_terms and len(classified_cases) == 0:
        logger.warning("No samples match the case phenotype terms.")
    if control_terms and len(classified_controls) == 0:
        logger.warning("No samples match the control phenotype terms.")

    logger.debug(f"Classified {len(classified_cases)} cases and {len(classified_controls)} controls from phenotypes.")
    return classified_cases, classified_controls

def build_sample_phenotype_map(df):
    """
    Build a sample->set_of_phenotypes map from the 'phenotypes' column if present.
    If 'phenotypes' is missing or empty, return {}.

    Logic:
    - For each variant line:
      * 'phenotypes' column may contain phenotype strings.
      * The phenotypes might be aligned with samples in the GT column if separated by ';'.
      * If multiple samples appear in GT, and phenotypes are also separated by ';',
        we align each phenotype subset with the corresponding sample.
    - We aggregate these phenotypes per sample across all variants.
    - This is a heuristic approach since phenotype data might not be fully variant-dependent.
    """

    if "phenotypes" not in df.columns:
        return {}

    sample_phenos = defaultdict(set)

    for idx, row in df.iterrows():
        pheno_str = row.get("phenotypes", "")
        if not isinstance(pheno_str, str) or not pheno_str.strip():
            continue

        # Extract the GT field and corresponding samples
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

        # Split phenotypes by ';' to see if we have multiple sets
        pheno_groups = pheno_str.split(";")
        pheno_groups = [pg.strip() for pg in pheno_groups if pg.strip()]

        if len(pheno_groups) > 1:
            # We have multiple phenotype groups, ideally one per sample
            if len(pheno_groups) == len(sample_names):
                # Perfect alignment: each phenotype group corresponds to one sample
                for sname, pgroup in zip(sample_names, pheno_groups):
                    phenos = {p.strip() for p in pgroup.split(",") if p.strip()}
                    sample_phenos[sname].update(phenos)
            else:
                # Mismatch in counts: log warning and fallback to assigning all phenotypes to all samples
                logger.warning(f"Number of phenotype groups ({len(pheno_groups)}) does not match sample count ({len(sample_names)}) at row {idx}. Assigning all phenotypes to all samples in this row.")
                # Combine all phenotypes from all groups
                combined_phenos = set()
                for pgroup in pheno_groups:
                    combined_phenos.update({p.strip() for p in pgroup.split(",") if p.strip()})
                for sname in sample_names:
                    sample_phenos[sname].update(combined_phenos)
        else:
            # Only one phenotype group (or none)
            # Assign all phenotypes to all samples in this row
            phenos = {p.strip() for p in pheno_str.split(",") if p.strip()}
            for sname in sample_names:
                sample_phenos[sname].update(phenos)

    return sample_phenos


def assign_case_control_counts(df, case_samples, control_samples, all_samples):
    """
    Create columns:
      - proband_count: total number of proband (case) samples (constant for all variants)
      - control_count: total number of control samples (constant for all variants)
      - proband_variant_count: how many proband samples carry the variant at this variant
      - control_variant_count: how many control samples carry the variant at this variant
      - proband_allele_count: total variant alleles among proband samples for this variant
      - control_allele_count: total variant alleles among control samples for this variant

    Logic:
    - proband_count/control_count are determined once from the sets case_samples/control_samples.
    - For each variant:
      * Parse GT column to find samples with variant alleles.
      * If a sample isn't listed, assume 0/0.
      * Count how many probands/controls have the variant and sum their allele counts.
    """

    logger.debug("Assigning case/control counts to variants...")

    # Determine total number of cases and controls once
    total_proband = len(case_samples)
    total_control = len(control_samples)
    logger.debug(f"Total proband samples: {total_proband}, Total control samples: {total_control}")

    proband_variant_count_list = []
    control_variant_count_list = []
    proband_allele_count_list = []
    control_allele_count_list = []

    # We'll store these totals as constant columns after the loop
    # proband_count and control_count remain constant per variant
    for idx, val in enumerate(df["GT"]):
        if idx % 1000 == 0:
            logger.debug(f"Processing variant row {idx} for allele counts...")

        # Map of samples in GT: sample_name -> genotype string
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

        # Iterate over all samples; those not in samples_with_variant are assumed 0/0
        for sample_name in all_samples:
            genotype = samples_with_variant.get(sample_name, "0/0")
            allele_count = genotype_to_allele_count(genotype)

            if sample_name in case_samples:
                # case sample
                if allele_count > 0:
                    p_variant_count += 1
                    p_allele_count += allele_count
            elif sample_name in control_samples:
                # control sample
                if allele_count > 0:
                    c_variant_count += 1
                    c_allele_count += allele_count
            # else excluded sample, ignore completely

        proband_variant_count_list.append(p_variant_count)
        control_variant_count_list.append(c_variant_count)
        proband_allele_count_list.append(p_allele_count)
        control_allele_count_list.append(c_allele_count)

    # Assign the columns
    # proband_count and control_count are constant for all variants:
    df["proband_count"] = total_proband
    df["control_count"] = total_control
    df["proband_variant_count"] = proband_variant_count_list
    df["control_variant_count"] = control_variant_count_list
    df["proband_allele_count"] = proband_allele_count_list
    df["control_allele_count"] = control_allele_count_list

    logger.debug("Case/control counts assigned. Example for first row:")
    if len(df) > 0:
        logger.debug(df.iloc[0][["proband_count","proband_variant_count","proband_allele_count","control_count","control_variant_count","control_allele_count"]].to_dict())

    return df


def extract_sample_and_genotype(sample_field):
    """
    Extract sample name and genotype from a field like 'sample(0/1)'.
    If parentheses are missing, assume no genotype is specified, treat as no variant.
    """
    start_idx = sample_field.find("(")
    end_idx = sample_field.find(")")
    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        sample_name = sample_field[:start_idx].strip()
        genotype = sample_field[start_idx+1:end_idx].strip()
        return sample_name, genotype
    else:
        # No genotype info provided, treat as empty genotype
        return sample_field.strip(), ""

def genotype_to_allele_count(genotype):
    """
    Convert genotype string to allele count:
    - '1/1' -> 2 alleles
    - '0/1' or '1/0' -> 1 allele
    - '0/0' or '' -> 0 alleles
    """
    if genotype == "1/1":
        return 2
    elif genotype in ["0/1", "1/0"]:
        return 1
    return 0

def compute_basic_stats(df, all_samples):
    """
    Compute basic statistics about the dataset:
    - Number of variants
    - Number of samples (from all_samples)
    - Number of genes
    - Het/Hom counts from GT fields
    - Variant type and impact counts if columns are present

    Debug logging included.
    """
    logger.debug("Computing basic stats...")
    num_variants = len(df)
    num_samples = len(all_samples)
    num_genes = df["GENE"].nunique()

    het_counts = 0
    hom_counts = 0
    for val in df["GT"]:
        if isinstance(val, str) and val.strip():
            for g in val.split(";"):
                g = g.strip()
                _, genotype = extract_sample_and_genotype(g)
                if genotype in ["0/1", "1/0"]:
                    het_counts += 1
                elif genotype == "1/1":
                    hom_counts += 1

    metric_rows = []
    metric_rows.append(["Number of variants", str(num_variants)])
    metric_rows.append(["Number of samples", str(num_samples)])
    metric_rows.append(["Number of genes", str(num_genes)])
    metric_rows.append(["Het counts", str(het_counts)])
    metric_rows.append(["Hom counts", str(hom_counts)])

    if "EFFECT" in df.columns:
        variant_types = df["EFFECT"].value_counts()
        for vt, count in variant_types.items():
            metric_rows.append([f"Variant_type_{vt}", str(count)])

    if "IMPACT" in df.columns:
        impact_types = df["IMPACT"].value_counts()
        for it, count in impact_types.items():
            metric_rows.append([f"Impact_type_{it}", str(count)])

    basic_stats_df = pd.DataFrame(metric_rows, columns=["metric", "value"])
    logger.debug("Basic stats computed.")
    return basic_stats_df

def compute_gene_stats(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute gene-level statistics by summing proband and control counts/alleles across variants.
    """
    logger.debug("Computing gene-level stats...")
    for col in ["proband_count", "control_count", "proband_allele_count", "control_allele_count"]:
        if col not in df.columns:
            df[col] = 0
    grouped = df.groupby("GENE").agg({
        "proband_count": "sum",
        "control_count": "sum",
        "proband_allele_count": "sum",
        "control_allele_count": "sum"
    }).reset_index()
    logger.debug("Gene-level stats computed.")
    return grouped

def compute_impact_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    If IMPACT column exists, compute per-gene impact summary.
    """
    logger.debug("Computing impact summary...")
    if "GENE" not in df.columns or "IMPACT" not in df.columns:
        logger.debug("No IMPACT or GENE column, returning empty impact summary.")
        return pd.DataFrame()
    impact_counts = df.groupby(["GENE", "IMPACT"]).size().reset_index(name="count")
    pivot_impact = impact_counts.pivot(index="GENE", columns="IMPACT", values="count").fillna(0).reset_index()
    logger.debug("Impact summary computed.")
    return pivot_impact

def compute_variant_type_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    If EFFECT column exists, compute per-gene variant type summary.
    """
    logger.debug("Computing variant type summary...")
    if "GENE" not in df.columns or "EFFECT" not in df.columns:
        logger.debug("No EFFECT or GENE column, returning empty variant type summary.")
        return pd.DataFrame()
    type_counts = df.groupby(["GENE", "EFFECT"]).size().reset_index(name="count")
    pivot_types = type_counts.pivot(index="GENE", columns="EFFECT", values="count").fillna(0).reset_index()
    logger.debug("Variant type summary computed.")
    return pivot_types

def merge_and_format_stats(gene_stats: pd.DataFrame,
                           impact_summary: pd.DataFrame,
                           variant_type_summary: pd.DataFrame) -> pd.DataFrame:
    """
    Merge gene_stats with impact_summary and variant_type_summary.
    """
    logger.debug("Merging gene stats with impact and variant type summaries...")
    merged = gene_stats.copy()
    if not impact_summary.empty:
        merged = pd.merge(merged, impact_summary, on="GENE", how="left")
    if not variant_type_summary.empty:
        merged = pd.merge(merged, variant_type_summary, on="GENE", how="left")
    merged = merged.fillna(0)
    logger.debug("All gene-level stats merged.")
    return merged

def gene_burden_fisher(subdf):
    """
    Perform Fisher's exact test for gene burden analysis on a subset of variants from a single gene.

    Columns used:
    - proband_allele_count
    - control_allele_count
    - proband_count
    - control_count

    Returns a DataFrame with fisher_p_value and related counts.
    """
    if "proband_allele_count" not in subdf.columns:
        subdf["proband_allele_count"] = 0
    if "control_allele_count" not in subdf.columns:
        subdf["control_allele_count"] = 0
    if "proband_count" not in subdf.columns:
        subdf["proband_count"] = 0
    if "control_count" not in subdf.columns:
        subdf["control_count"] = 0

    proband_alleles = subdf["proband_allele_count"].sum()
    control_alleles = subdf["control_allele_count"].sum()
    max_proband_count = subdf["proband_count"].max()
    max_control_count = subdf["control_count"].max()

    total_case_samples = max_proband_count if max_proband_count > 0 else 0
    total_control_samples = max_control_count if max_control_count > 0 else 0

    if total_case_samples == 0 and total_control_samples == 0:
        return pd.DataFrame([{
            "GENE": subdf["GENE"].iloc[0],
            "proband_alleles": 0,
            "control_alleles": 0,
            "max_proband_count": 0,
            "max_control_count": 0,
            "proband_ref_alleles": 0,
            "control_ref_alleles": 0,
            "fisher_p_value": 1.0
        }])

    proband_ref_alleles = (total_case_samples * 2) - proband_alleles if total_case_samples else 0
    control_ref_alleles = (total_control_samples * 2) - control_alleles if total_control_samples else 0

    table = [[proband_alleles, control_alleles],
             [proband_ref_alleles, control_ref_alleles]]

    # Log the 2x2 table for debugging
    logger.debug(f"Fisher test table for gene {subdf['GENE'].iloc[0]}: {table}")

    oddsratio, p_value = fisher_exact(table) if fisher_exact else (None, 1.0)

    return pd.DataFrame([{
        "GENE": subdf["GENE"].iloc[0],
        "proband_alleles": proband_alleles,
        "control_alleles": control_alleles,
        "max_proband_count": max_proband_count,
        "max_control_count": max_control_count,
        "proband_ref_alleles": proband_ref_alleles,
        "control_ref_alleles": control_ref_alleles,
        "fisher_p_value": p_value
    }])

# File: variantcentrifuge/analyze_variants.py
# Location: variantcentrifuge/variantcentrifuge/analyze_variants.py

"""
Variant analysis module for gene burden and other statistics.

This module:
- Reads a TSV of variants and their annotations (including a GT column listing sample genotypes).
- Classifies samples into case/control sets based on user input (sample lists or phenotype terms).
- If no case/control criteria are provided, defaults to making all samples controls.
- Computes variant-level statistics and, if requested, performs a gene burden analysis (per-gene) via Fisher's exact test.

Fallback Behavior:
- If case or control groups are not provided by any method (phenotypes, sample lists), 
  all samples are considered controls by default.

Logging:
- Logs at each major step and data processing point.
- Logs sizes of sample sets, classification outcomes, and stats computations.

Comments:
- Each function has a docstring explaining logic and parameters.
- Inline comments clarify specific code logic.
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

import statsmodels.stats.multitest as smm

logger = logging.getLogger("variantcentrifuge")


def analyze_variants(lines, cfg):
    """
    Analyze variants and optionally perform gene burden analysis.

    Steps:
    - Parse input TSV into a DataFrame.
    - Retrieve full sample list from cfg["sample_list"].
    - Determine case/control sets based on cfg (samples or phenotypes).
    - Compute per-variant case/control allele counts.
    - Compute basic and optionally comprehensive gene-level stats.
    - If requested, perform gene burden (Fisher's exact test + correction).

    Parameters
    ----------
    lines : iterator of str
        Input lines representing a TSV with variant data.
    cfg : dict
        Configuration dictionary with keys:
        - sample_list (str): Comma-separated full sample list from VCF
        - case_samples, control_samples (list of str)
        - case_phenotypes, control_phenotypes (list of str)
        - perform_gene_burden (bool): Whether to perform gene burden analysis
        - gene_burden_mode (str): "samples" or "alleles"
        - correction_method (str): "fdr" or "bonferroni"
        - no_stats (bool): Skip stats if True
        - stats_output_file (str): Path to stats output
        - gene_burden_output_file (str, optional)
        - xlsx (bool): If True, might append to Excel after analysis

    Yields
    ------
    str
        Processed lines of output TSV or gene-level burden results.
    """
    perform_gene_burden = cfg.get("perform_gene_burden", False)
    stats_output_file = cfg.get("stats_output_file")
    no_stats = cfg.get("no_stats", False)

    logger.debug(f"analyze_variants: perform_gene_burden={perform_gene_burden}, "
                 f"stats_output_file={stats_output_file}, no_stats={no_stats}")

    # Read all input lines into a single text block
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

    # Ensure we have the full sample list from cfg
    if "sample_list" not in cfg or not cfg["sample_list"].strip():
        logger.error("No sample_list found in cfg. Unable to determine the full sample set.")
        for line in text_data.strip().split("\n"):
            yield line
        return

    # all_samples: full set of samples
    all_samples = set(cfg["sample_list"].split(","))
    logger.debug(f"Total samples from VCF header: {len(all_samples)}")

    # Determine case/control sets
    case_samples, control_samples = determine_case_control_sets(all_samples, cfg, df)
    logger.debug(f"Number of case samples: {len(case_samples)}; Number of control samples: {len(control_samples)}")

    # Assign case/control counts per variant
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

    # Perform gene burden analysis if requested
    if perform_gene_burden:
        logger.debug("Performing configurable gene burden analysis...")
        gene_burden_results = perform_gene_burden_analysis(df, cfg)

        if gene_burden_results.empty:
            logger.warning("Gene burden results are empty. Returning variant-level data.")
            out_str = df.to_csv(sep="\t", index=False)
            for line in out_str.strip().split("\n"):
                yield line
            return

        # Write burden results to file if provided
        burden_output_file = cfg.get("gene_burden_output_file", "gene_burden_results.tsv")
        logger.info(f"Writing gene burden results to {burden_output_file}")
        gene_burden_results.to_csv(burden_output_file, sep="\t", index=False)

        # Yield burden results as output
        burden_text = gene_burden_results.to_csv(sep="\t", index=False)
        for line in burden_text.strip().split("\n"):
            yield line
    else:
        logger.debug("Gene burden not requested, returning processed data.")
        out_str = df.to_csv(sep="\t", index=False)
        for line in out_str.strip().split("\n"):
            yield line


def perform_gene_burden_analysis(df, cfg):
    """
    Perform gene burden analysis for each gene.

    Steps:
    1. Aggregate variant counts per gene based on chosen mode (samples or alleles).
    2. Perform Fisher's exact test for each gene.
    3. Apply multiple testing correction (FDR or Bonferroni).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with per-variant case/control counts.
    cfg : dict
        Configuration dictionary with keys:
        - gene_burden_mode: 'samples' or 'alleles'
        - correction_method: 'fdr' or 'bonferroni'

    Returns
    -------
    pd.DataFrame
        A DataFrame with gene-level burden results including p-values and odds ratios.
    """
    logger.debug("Starting gene burden aggregation...")

    mode = cfg.get("gene_burden_mode", "alleles")
    correction_method = cfg.get("correction_method", "fdr")

    # Aggregate per-gene counts
    grouped = df.groupby("GENE").agg({
        "proband_count": "max",
        "control_count": "max",
        "proband_variant_count": "sum",
        "control_variant_count": "sum",
        "proband_allele_count": "sum",
        "control_allele_count": "sum"
    }).reset_index()

    results = []
    for _, row in grouped.iterrows():
        gene = row["GENE"]
        p_count = int(row["proband_count"])
        c_count = int(row["control_count"])

        # Skip if no samples in either group
        if p_count == 0 and c_count == 0:
            continue

        if mode == "samples":
            # Per-sample variant counts
            p_var = int(row["proband_variant_count"])
            c_var = int(row["control_variant_count"])
            p_ref = p_count - p_var
            c_ref = c_count - c_var
            table = [[p_var, c_var],
                     [p_ref, c_ref]]
            var_metric = ("proband_variant_count", "control_variant_count")
        else:
            # Per-allele counts
            p_all = int(row["proband_allele_count"])
            c_all = int(row["control_allele_count"])
            p_ref = p_count * 2 - p_all
            c_ref = c_count * 2 - c_all
            table = [[p_all, c_all],
                     [p_ref, c_ref]]
            var_metric = ("proband_allele_count", "control_allele_count")

        # Fisher's exact test if available
        if fisher_exact is not None:
            odds_ratio, pval = fisher_exact(table)
        else:
            odds_ratio = None
            pval = 1.0

        results.append({
            "GENE": gene,
            "proband_count": p_count,
            "control_count": c_count,
            var_metric[0]: table[0][0],
            var_metric[1]: table[0][1],
            "raw_p_value": pval,
            "odds_ratio": odds_ratio
        })

    if not results:
        logger.warning("No genes found with variant data for gene burden analysis.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)

    # Multiple testing correction
    pvals = results_df["raw_p_value"].values
    if correction_method == "bonferroni":
        corrected_pvals = smm.multipletests(pvals, method="bonferroni")[1]
    else:
        # Default to FDR (Benjamini–Hochberg)
        corrected_pvals = smm.multipletests(pvals, method="fdr_bh")[1]

    results_df["corrected_p_value"] = corrected_pvals

    if mode == "samples":
        final_cols = [
            "GENE",
            "proband_count", "control_count",
            "proband_variant_count", "control_variant_count",
            "raw_p_value", "corrected_p_value", "odds_ratio"
        ]
    else:
        final_cols = [
            "GENE",
            "proband_count", "control_count",
            "proband_allele_count", "control_allele_count",
            "raw_p_value", "corrected_p_value", "odds_ratio"
        ]

    results_df = results_df[final_cols]
    logger.debug("Gene burden analysis complete.")
    return results_df


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
            # Only case samples given
            control_samples = all_samples - case_samples
        elif control_samples and not case_samples:
            # Only control samples given
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


def compute_basic_stats(df, all_samples):
    """
    Compute basic statistics about the dataset:
    - Number of variants
    - Number of samples
    - Number of genes
    - Het/Hom counts from GT fields
    - Variant type and impact counts if available

    Parameters
    ----------
    df : pd.DataFrame
        Input variants DataFrame.
    all_samples : set
        Set of all samples.

    Returns
    -------
    pd.DataFrame
        A DataFrame with 'metric' and 'value' columns listing basic stats.
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

    metric_rows = [
        ["Number of variants", str(num_variants)],
        ["Number of samples", str(num_samples)],
        ["Number of genes", str(num_genes)],
        ["Het counts", str(het_counts)],
        ["Hom counts", str(hom_counts)]
    ]

    # Count variant and impact types if columns present
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
    Compute gene-level aggregated stats (sum of proband/control counts and alleles).

    Parameters
    ----------
    df : pd.DataFrame
        Input variants DataFrame with assigned case/control counts.

    Returns
    -------
    pd.DataFrame
        Gene-level summary DataFrame.
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
    Compute per-gene impact summary if IMPACT column exists.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame

    Returns
    -------
    pd.DataFrame
        A pivoted table of gene vs impact counts.
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
    Compute per-gene variant type summary if EFFECT column exists.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame

    Returns
    -------
    pd.DataFrame
        A pivoted table of gene vs variant types.
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
    Merge gene_stats with impact_summary and variant_type_summary into a single DataFrame.

    Parameters
    ----------
    gene_stats : pd.DataFrame
        DataFrame of gene-level aggregated stats
    impact_summary : pd.DataFrame
        DataFrame of gene vs impact counts
    variant_type_summary : pd.DataFrame
        DataFrame of gene vs variant type counts

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with all gene-level stats.
    """
    logger.debug("Merging gene stats with impact and variant type summaries...")
    merged = gene_stats.copy()
    if not impact_summary.empty:
        merged = pd.merge(merged, impact_summary, on="GENE", how="left")
    if not variant_type_summary.empty:
        merged = pd.merge(merged, variant_type_summary, on="GENE", how="left")
    # Fill NaN with 0 for counts
    merged = merged.fillna(0)
    logger.debug("All gene-level stats merged.")
    return merged

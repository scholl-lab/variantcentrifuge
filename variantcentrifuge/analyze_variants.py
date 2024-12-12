# File: variantcentrifuge/analyze_variants.py
# Location: variantcentrifuge/variantcentrifuge/analyze_variants.py

"""
Variant analysis module.

This module performs basic variant statistics and an optional gene burden analysis
using Fisher's exact test, replicating the logic from the analyze_variants.R script.
If gene burden analysis is requested, it returns the burden results instead of the
original lines. Otherwise, it can pass through the original lines.

It can also write statistics to a separate file if configured.
"""

import io
import pandas as pd
from collections import Counter
from .utils import log_message
from math import isnan

try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

def analyze_variants(lines, cfg):
    """
    Analyze variants, potentially performing gene burden analysis and writing statistics.

    Parameters
    ----------
    lines : iterator
        Iterator of TSV lines containing variant data.
    cfg : dict
        Configuration dictionary, which may contain:
        - perform_gene_burden (bool): Whether to perform gene burden analysis.
        - stats_output_file (str): Path to write statistics file.
        - Additional assumptions: The input lines contain at least the required columns.

    Returns
    -------
    iterator
        If perform_gene_burden is True, yields the gene burden analysis results as TSV lines.
        Otherwise, yields the original lines unchanged.
    """

    perform_gene_burden = cfg.get("perform_gene_burden", False)
    stats_output_file = cfg.get("stats_output_file")

    # Read all lines into a DataFrame
    # We assume lines are TSV with a header.
    text_data = "".join(line for line in lines)
    if not text_data.strip():
        # No data
        return

    df = pd.read_csv(io.StringIO(text_data), sep="\t", dtype=str)

    # Convert numeric columns
    numeric_cols = ["proband_count", "proband_allele_count", "control_count", "control_allele_count"]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
        else:
            log_message("ERROR", f"Missing required column: {col}")
            return

    required_columns = ["CHROM", "POS", "REF", "ALT", "GENE", "GT",
                        "proband_count", "proband_allele_count",
                        "control_count", "control_allele_count"]
    missing_columns = [c for c in required_columns if c not in df.columns]
    if missing_columns:
        log_message("ERROR", f"Missing required columns: {', '.join(missing_columns)}")
        return

    # Calculate basic statistics
    # Number of variants
    num_variants = len(df)

    # Number of samples: unique samples from GT field
    # GT field contains samples separated by ";" or empty if no genotype?
    # According to replacer logic, GT might contain multiple samples separated by ";"
    # We'll split all non-empty GT fields and get unique samples.
    all_samples = set()
    for val in df["GT"]:
        if isinstance(val, str) and val.strip():
            for s in val.split(";"):
                s = s.strip()
                if s:
                    # If sample appended genotype e.g. "sample(0/1)" remove that part if needed
                    # The replacer script may have appended genotypes in parentheses.
                    # We'll strip that off.
                    # sample(0/1) => sample
                    idx = s.find("(")
                    if idx != -1:
                        s = s[:idx]
                    all_samples.add(s)
    num_samples = len(all_samples)

    # Number of genes
    num_genes = df["GENE"].nunique()

    # Het and Hom counts:
    # Hets: Count occurrences of "0/1" or "1/0" in GT field
    # Homs: Count occurrences of "1/1"
    # Note: GT might contain multiple samples separated by ";"
    # We must count occurrences across all samples in that field.
    het_counts = 0
    hom_counts = 0
    for val in df["GT"]:
        if isinstance(val, str) and val.strip():
            for g in val.split(";"):
                g = g.strip()
                # Extract genotype if appended
                # sample(0/1) => 0/1
                genotype = extract_genotype(g)
                if genotype == "0/1" or genotype == "1/0":
                    het_counts += 1
                elif genotype == "1/1":
                    hom_counts += 1

    # Variant types by EFFECT (if EFFECT column exists)
    variant_types = None
    if "EFFECT" in df.columns:
        variant_types = df["EFFECT"].value_counts().reset_index()
        variant_types.columns = ["EFFECT", "count"]

    # Impact types by IMPACT (if IMPACT column exists)
    impact_types = None
    if "IMPACT" in df.columns:
        impact_types = df["IMPACT"].value_counts().reset_index()
        impact_types.columns = ["IMPACT", "count"]

    # Print/log some stats
    log_message("INFO", f"Number of variants: {num_variants}")
    log_message("INFO", f"Number of samples: {num_samples}")
    log_message("INFO", f"Number of genes: {num_genes}")
    log_message("INFO", f"Het counts: {het_counts}")
    log_message("INFO", f"Hom counts: {hom_counts}")

    if variant_types is not None:
        log_message("INFO", f"Variant types:\n{variant_types}")
    else:
        log_message("INFO", "No EFFECT column present, skipping variant types.")

    if impact_types is not None:
        log_message("INFO", f"Impact types:\n{impact_types}")
    else:
        log_message("INFO", "No IMPACT column present, skipping impact types.")

    # If perform_gene_burden is True, do gene burden analysis
    if perform_gene_burden:
        if fisher_exact is None:
            log_message("ERROR", "scipy not available for Fisher test.")
            return

        # perform gene burden analysis
        # Group by GENE and calculate:
        # proband_alleles = sum(proband_allele_count)
        # control_alleles = sum(control_allele_count)
        # max_proband_count = max(proband_count)
        # max_control_count = max(control_count)
        # proband_ref_alleles = sum(max_proband_count*2 - proband_allele_count)
        # control_ref_alleles = sum(max_control_count*2 - control_allele_count)
        # fisher_p_value = fisher.test(...)
        grouped = df.groupby("GENE").apply(gene_burden_fisher)
        # grouped is now a DataFrame with columns defined in gene_burden_fisher function

        burden_text = grouped.to_csv(sep="\t", index=False)

        # If stats file is specified, write stats now
        if stats_output_file:
            write_stats(stats_output_file, num_variants, num_samples, num_genes,
                        het_counts, hom_counts, variant_types, impact_types)

        # yield burden analysis lines
        for line in burden_text.strip().split("\n"):
            yield line
    else:
        # If not performing gene burden, just return original lines
        # If stats_output_file is given, write stats
        if stats_output_file:
            write_stats(stats_output_file, num_variants, num_samples, num_genes,
                        het_counts, hom_counts, variant_types, impact_types)

        # Return original data
        for line in text_data.strip().split("\n"):
            yield line


def extract_genotype(sample_field):
    """
    Extract the genotype from a field like 'sample(0/1)' or 'sample'.
    If no parentheses, assume no genotype appended.
    """
    # If we have parentheses
    start_idx = sample_field.find("(")
    end_idx = sample_field.find(")")
    if start_idx != -1 and end_idx != -1 and end_idx > start_idx:
        return sample_field[start_idx+1:end_idx]
    # If no parentheses, might be just a sample name or empty
    # No genotype info if missing parentheses. In original logic,
    # if no genotype appended and not "0/0", we might not know genotype?
    # We'll assume if no parentheses, no genotype info. Return empty or something.
    # Or assume original replacer ensures genotype appended if needed.
    # To be safe, return empty string.
    return ""


def gene_burden_fisher(subdf):
    """
    Perform gene burden calculation for a group of variants from one gene.

    subdf: DataFrame filtered for a single gene.
    Returns: A single-row DataFrame with columns:
     GENE, proband_alleles, control_alleles,
     max_proband_count, max_control_count,
     proband_ref_alleles, control_ref_alleles,
     fisher_p_value
    """
    proband_alleles = subdf["proband_allele_count"].sum()
    control_alleles = subdf["control_allele_count"].sum()
    max_proband_count = subdf["proband_count"].max()
    max_control_count = subdf["control_count"].max()

    # compute ref alleles
    proband_ref_alleles = (max_proband_count * 2 * len(subdf)) - proband_alleles
    control_ref_alleles = (max_control_count * 2 * len(subdf)) - control_alleles

    # fisher test
    table = [[proband_alleles, control_alleles],
             [proband_ref_alleles, control_ref_alleles]]

    oddsratio, p_value = fisher_exact(table)

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


def write_stats(stats_file, num_variants, num_samples, num_genes, het_counts, hom_counts,
                variant_types, impact_types):
    """
    Write statistics to a file.
    """
    with open(stats_file, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"Number of variants\t{num_variants}\n")
        f.write(f"Number of samples\t{num_samples}\n")
        f.write(f"Number of genes\t{num_genes}\n")
        f.write(f"Het counts\t{het_counts}\n")
        f.write(f"Hom counts\t{hom_counts}\n")

        if variant_types is not None:
            f.write("Variant types:\n")
            for _, row in variant_types.iterrows():
                f.write(f"{row['EFFECT']}\t{row['count']}\n")

        if impact_types is not None:
            f.write("Impact types:\n")
            for _, row in impact_types.iterrows():
                f.write(f"{row['IMPACT']}\t{row['count']}\n")

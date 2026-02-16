# File: variantcentrifuge/phenotype.py
# Location: variantcentrifuge/variantcentrifuge/phenotype.py

"""
Phenotype integration module.

This module loads a phenotype file containing sample-to-phenotype mappings and provides a function
to aggregate phenotypes for a given list of samples.

- The phenotype file must be .csv or .tsv (detected by extension).
- The specified sample and phenotype columns must be present in the file.
- Phenotypes are stored in a dictionary (sample -> set of phenotypes).
- Given a list of samples, phenotypes are aggregated as follows:
  - For each sample, join multiple phenotypes by ",".
  - For multiple samples, join each sample's phenotype string by ";".
"""

import logging
import os

logger = logging.getLogger("variantcentrifuge")


def load_phenotypes(
    phenotype_file: str, sample_column: str, phenotype_column: str
) -> dict[str, set[str]]:
    """
    Load phenotypes from a .csv or .tsv file into a dictionary.

    Parameters
    ----------
    phenotype_file : str
        Path to the phenotype file (must be .csv or .tsv).
    sample_column : str
        Name of the column containing sample IDs.
    phenotype_column : str
        Name of the column containing phenotype values.

    Returns
    -------
    dict of {str: set of str}
        A dictionary mapping each sample to a set of associated phenotypes.

    Raises
    ------
    ValueError
        If the file is not .csv or .tsv, or if the required columns are not found.
    """
    if phenotype_file.endswith(".csv"):
        delim = ","
    elif phenotype_file.endswith(".tsv"):
        delim = "\t"
    else:
        raise ValueError("Phenotype file must be .csv or .tsv")

    phenotypes: dict[str, set[str]] = {}
    with open(phenotype_file, encoding="utf-8") as f:
        header_line = f.readline().rstrip("\n")
        header = header_line.split(delim)

        if sample_column not in header:
            raise ValueError(f"Sample column '{sample_column}' not found in phenotype file.")
        if phenotype_column not in header:
            raise ValueError(f"Phenotype column '{phenotype_column}' not found in phenotype file.")

        sample_idx = header.index(sample_column)
        pheno_idx = header.index(phenotype_column)

        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split(delim)
            if len(fields) <= max(sample_idx, pheno_idx):
                continue
            samp = fields[sample_idx].strip()
            pheno = fields[pheno_idx].strip()
            if samp not in phenotypes:
                phenotypes[samp] = set()
            if pheno:
                phenotypes[samp].add(pheno)

    return phenotypes


def aggregate_phenotypes_for_samples(samples: list[str], phenotypes: dict[str, set[str]]) -> str:
    """
    Aggregate phenotypes for a given list of samples into a single string.

    For each sample:
    - Join multiple phenotypes with ",".
    For multiple samples:
    - Join each sample's phenotype string with ";".

    Parameters
    ----------
    samples : list of str
        List of sample IDs.
    phenotypes : dict of {str: set of str}
        Dictionary mapping sample IDs to a set of phenotypes.

    Returns
    -------
    str
        A string aggregating all phenotypes for the given samples, with phenotypes
        comma-separated per sample, and samples separated by ";".
    """
    sample_phenos = []
    for s in samples:
        p_str = ",".join(sorted(phenotypes[s])) if phenotypes.get(s) else ""
        sample_phenos.append(p_str)

    return ";".join(sample_phenos)


def format_phenotypes_like_gt_column(samples: list[str], phenotypes: dict[str, set[str]]) -> str:
    """
    Format phenotypes in the same style as the GT column with sample IDs.

    Creates a string similar to GT column format:
    "SampleID(phenotype1,phenotype2);SampleID(phenotype3);..."

    This matches the format used in genotype replacement where each sample's data is
    prefixed with the sample ID in parentheses.

    Parameters
    ----------
    samples : list of str
        List of sample IDs in the order they appear in the VCF.
    phenotypes : dict of {str: set of str}
        Dictionary mapping sample IDs to a set of phenotypes.

    Returns
    -------
    str
        A string with phenotypes formatted like GT column:
        "Sample1(pheno1,pheno2);Sample2(pheno3);..."
        Samples without phenotypes get empty parentheses: "Sample3()".
    """
    sample_entries = []
    for sample_id in samples:
        phenotype_str = ",".join(sorted(phenotypes[sample_id])) if phenotypes.get(sample_id) else ""

        # Format like GT column: SampleID(phenotypes)
        sample_entries.append(f"{sample_id}({phenotype_str})")

    return ";".join(sample_entries)


def extract_phenotypes_for_gt_row(gt_value: str, phenotypes: dict[str, set[str]]) -> str:
    """
    Extract phenotypes for samples that have variants in a specific GT row.

    Parses the GT column value to find which samples have variants, then returns
    their phenotypes in the same format as the GT column.

    Parameters
    ----------
    gt_value : str
        GT column value like "Sample1(0/1);Sample2(1/1);Sample3(./.)"
    phenotypes : dict of {str: set of str}
        Dictionary mapping sample IDs to a set of phenotypes.

    Returns
    -------
    str
        Phenotypes for samples with variants: "Sample1(pheno1,pheno2);Sample2(pheno3)"
        Samples with no variants (./. or 0/0) are excluded.
    """
    if not gt_value or not gt_value.strip():
        return ""

    phenotype_entries = []

    # Parse GT column to extract samples with variants
    for sample_entry in gt_value.split(";"):
        sample_entry = sample_entry.strip()
        if not sample_entry:
            continue

        # Extract sample ID and genotype: "Sample1(0/1)" -> "Sample1", "0/1"
        if "(" in sample_entry and sample_entry.endswith(")"):
            sample_id = sample_entry.split("(")[0]
            genotype = sample_entry.split("(")[1][:-1]  # Remove closing )

            # Skip entries without sample ID (malformed entries like "(0/1)")
            if not sample_id or not sample_id.strip():
                continue

            # Skip samples with no variant (./. or 0/0)
            if genotype in ["./.", "0/0", ""]:
                continue

            # Get phenotypes for this sample
            if phenotypes.get(sample_id):
                phenotype_str = ",".join(sorted(phenotypes[sample_id]))
                phenotype_entries.append(f"{sample_id}({phenotype_str})")
            else:
                # Sample has variant but no phenotypes
                phenotype_entries.append(f"{sample_id}()")

    return ";".join(phenotype_entries)


def extract_phenotypes_from_sample_columns(
    row, vcf_samples: list[str], phenotypes: dict[str, set[str]]
) -> str:
    """
    Extract phenotypes from per-sample GT columns (Phase 11 bcftools output).

    With bcftools query output, each sample has its own GT column (GEN[0].GT, GEN[1].GT, ...).
    This function examines these columns to find samples with variants and returns their phenotypes.

    Parameters
    ----------
    row : pd.Series or namedtuple
        DataFrame row containing per-sample GT columns
    vcf_samples : list of str
        List of sample IDs in VCF order (matches GEN[N] indexing)
    phenotypes : dict of {str: set of str}
        Dictionary mapping sample IDs to a set of phenotypes.

    Returns
    -------
    str
        Phenotypes for samples with variants: "Sample1(pheno1,pheno2);Sample2(pheno3)"
        Samples with reference genotypes (0/0, ./., NA, empty) are excluded.
    """
    import pandas as pd

    phenotype_entries = []

    for sample_id in vcf_samples:
        # Get GT value for this sample - try both sanitized and original column names
        # Original: "GEN[0].GT", Sanitized: "GEN_0__GT"
        gt_value = None

        # Try getting by attribute (namedtuple or Series)
        if hasattr(row, sample_id):
            gt_value = getattr(row, sample_id)
        # Try dictionary-style access (Series)
        elif hasattr(row, "__getitem__"):
            import contextlib

            with contextlib.suppress(KeyError, IndexError):
                gt_value = row[sample_id]

        # Handle missing/NA/reference genotypes
        if pd.isna(gt_value) or not gt_value:
            continue

        gt_str = str(gt_value).strip()
        if not gt_str or gt_str in ["0/0", "./.", ".", "NA"]:
            continue

        # Sample has variant - get phenotypes
        if phenotypes.get(sample_id):
            phenotype_str = ",".join(sorted(phenotypes[sample_id]))
            phenotype_entries.append(f"{sample_id}({phenotype_str})")
        else:
            # Sample has variant but no phenotypes
            phenotype_entries.append(f"{sample_id}()")

    return ";".join(phenotype_entries)

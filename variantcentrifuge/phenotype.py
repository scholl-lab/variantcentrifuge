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

import os
import logging
from typing import Dict, Set, List

logger = logging.getLogger("variantcentrifuge")


def load_phenotypes(phenotype_file: str, sample_column: str, phenotype_column: str) -> Dict[str, Set[str]]:
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

    phenotypes: Dict[str, Set[str]] = {}
    with open(phenotype_file, "r", encoding="utf-8") as f:
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


def aggregate_phenotypes_for_samples(samples: List[str], phenotypes: Dict[str, Set[str]]) -> str:
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
        if s in phenotypes and phenotypes[s]:
            p_str = ",".join(sorted(phenotypes[s]))
        else:
            p_str = ""
        sample_phenos.append(p_str)

    return ";".join(sample_phenos)

# File: variantcentrifuge/phenotype.py
# Location: variantcentrifuge/variantcentrifuge/phenotype.py

"""
Phenotype integration module.

This module loads a phenotype file containing sample-to-phenotype mappings and provides a function
to aggregate phenotypes for a given list of samples.

- The phenotype file must be .csv or .tsv (detected by extension).
- We find the specified sample column and phenotype column by name.
- Store phenotypes in a dictionary: sample -> set of phenotypes.
- Given a list of samples, aggregate their phenotypes into a single string.
  - For each sample, join multiple phenotypes by ",".
  - Join multiple samples by ";" in the final output.
"""

import os
import logging

logger = logging.getLogger("variantcentrifuge")

def load_phenotypes(phenotype_file, sample_column, phenotype_column):
    # Determine delimiter
    if phenotype_file.endswith(".csv"):
        delim = ","
    elif phenotype_file.endswith(".tsv"):
        delim = "\t"
    else:
        raise ValueError("Phenotype file must be .csv or .tsv")

    phenotypes = {}
    with open(phenotype_file, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split(delim)
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

def aggregate_phenotypes_for_samples(samples, phenotypes):
    """
    Given a list of samples and a phenotypes dict: sample -> set of phenotypes,
    return a single string representing all phenotypes for these samples.

    - For each sample, join multiple phenotypes with ",".
    - Then join multiple samples with ";".
    - If a sample has no phenotypes, produce empty string for that sample.
    """
    sample_phenos = []
    for s in samples:
        if s in phenotypes and phenotypes[s]:
            p_str = ",".join(sorted(phenotypes[s]))
        else:
            p_str = ""
        sample_phenos.append(p_str)

    # Join multiple samples with ";"
    return ";".join(sample_phenos)

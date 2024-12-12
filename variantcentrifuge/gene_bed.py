# File: variantcentrifuge/gene_bed.py
# Location: variantcentrifuge/variantcentrifuge/gene_bed.py

"""
Gene BED extraction module.

This module provides functionality to run snpEff genes2bed and produce
a BED file representing genes of interest. It also handles sorting and
optionally adding a 'chr' prefix to chromosome names.
"""

import subprocess
import tempfile
import os
from .utils import log_message


def get_gene_bed(reference, gene_name, interval_expand=0, add_chr=True):
    """
    Generate a BED file for the given gene(s) using snpEff genes2bed.
    Then sort and optionally add 'chr' prefix to chromosome names.

    Parameters
    ----------
    reference : str
        The reference name compatible with snpEff.
    gene_name : str
        Name of the gene or 'all' to consider all genes.
    interval_expand : int, optional
        Number of bases to expand the gene interval upstream/downstream.
    add_chr : bool, optional
        Whether to add 'chr' prefix to chromosome names.

    Returns
    -------
    str
        Path to the final BED file.
    """
    log_message("DEBUG", f"Entering get_gene_bed with reference={reference}, "
                         f"gene_name={gene_name}, interval_expand={interval_expand}, "
                         f"add_chr={add_chr}")

    # Create a temporary BED file
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "genes2bed", reference]
    if gene_name and gene_name.lower() != "all":
        cmd.append(gene_name)
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])

    log_message("INFO", f"Running: {' '.join(cmd)}")
    log_message("DEBUG", "Executing snpEff genes2bed command...")
    # Run snpEff genes2bed and capture output in bed_path
    subprocess.run(cmd, stdout=open(bed_path, "w", encoding="utf-8"), check=True)
    log_message("DEBUG", f"snpEff genes2bed completed, BED file at {bed_path}")

    # Sort the BED file
    sorted_bed = bed_path + ".sorted"
    sort_cmd = ["sortBed", "-i", bed_path]
    log_message("DEBUG", f"Sorting BED file with: {' '.join(sort_cmd)}")
    subprocess.run(sort_cmd, stdout=open(sorted_bed, "w", encoding="utf-8"), check=True)
    log_message("DEBUG", f"BED sorting completed, sorted file at {sorted_bed}")

    # Optionally add 'chr' prefix
    if add_chr:
        chr_bed = sorted_bed + ".chr"
        log_message("DEBUG", f"Adding 'chr' prefix to BED file {sorted_bed}")
        with open(chr_bed, "w", encoding="utf-8") as out_f, \
             open(sorted_bed, "r", encoding="utf-8") as in_f:
            for line in in_f:
                if not line.startswith("chr"):
                    out_f.write("chr" + line)
                else:
                    out_f.write(line)
        os.remove(sorted_bed)
        log_message("DEBUG", f"Final BED file with 'chr': {chr_bed}")
        return chr_bed

    log_message("DEBUG", f"Final BED file: {sorted_bed}")
    return sorted_bed

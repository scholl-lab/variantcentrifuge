# File: variantcentrifuge/gene_bed.py
# Location: variantcentrifuge/variantcentrifuge/gene_bed.py

"""
Gene BED extraction module.

This module provides functionality to run snpEff genes2bed and produce
a BED file representing genes of interest. It also handles sorting and
optional prefixing of chromosome names.
"""

import subprocess
import tempfile
import os
from .utils import log_message


def get_gene_bed(reference, gene_name, interval_expand=0, add_chr=True):
    """
    Generate a BED file for the given gene(s) using snpEff genes2bed.

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
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "genes2bed", reference]
    if gene_name and gene_name.lower() != "all":
        cmd.append(gene_name)
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])

    log_message("INFO", f"Running: {' '.join(cmd)}")
    with open(bed_path, "w", encoding="utf-8") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

    sorted_bed = bed_path + ".sorted"
    sort_cmd = ["sortBed", "-i", bed_path]
    with open(sorted_bed, "w", encoding="utf-8") as out_f:
        subprocess.run(sort_cmd, stdout=out_f, check=True)

    if add_chr:
        chr_bed = sorted_bed + ".chr"
        with open(chr_bed, "w", encoding="utf-8") as out_f, \
             open(sorted_bed, "r", encoding="utf-8") as in_f:
            for line in in_f:
                if not line.startswith("chr"):
                    out_f.write("chr" + line)
                else:
                    out_f.write(line)
        os.remove(sorted_bed)
        return chr_bed

    return sorted_bed

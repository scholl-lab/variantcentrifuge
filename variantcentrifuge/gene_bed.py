# File: variantcentrifuge/gene_bed.py
# Location: variantcentrifuge/variantcentrifuge/gene_bed.py

"""
Gene BED extraction module.

This module provides functionality to run snpEff genes2bed and produce
a BED file representing genes of interest. It also handles sorting and
optionally adding a 'chr' prefix to chromosome names.

Caching logic:
- If genes == "all", run snpEff genes2bed without specifying any genes or file.
- If multiple genes are provided, write them to a file (one gene per line) in the cache directory
  and pass the file to snpEff genes2bed using '-f <file>'.
- Cache results based on a hash of parameters and gene list.
"""

import subprocess
import tempfile
import os
import hashlib
import shutil
from .utils import log_message

def get_gene_bed(reference, gene_name, interval_expand=0, add_chr=True, output_dir="output"):
    """
    Generate a BED file for the given gene(s) using snpEff genes2bed.
    If gene_name == "all", run once with no gene arguments.
    If multiple genes are provided, write them to a file and use snpEff genes2bed -f <file>.

    Parameters
    ----------
    reference : str
        The reference name compatible with snpEff.
    gene_name : str
        "all" or space-separated list of genes.
    interval_expand : int, optional
        Number of bases to expand the gene interval upstream/downstream.
    add_chr : bool, optional
        Whether to add 'chr' prefix to chromosome names.
    output_dir : str, optional
        Directory where cached BED files can be stored.

    Returns
    -------
    str
        Path to the final BED file.
    """
    log_message("DEBUG", f"Entering get_gene_bed with reference={reference}, "
                         f"gene_name={gene_name}, interval_expand={interval_expand}, "
                         f"add_chr={add_chr}")

    cache_dir = os.path.join(output_dir, "bed_cache")
    os.makedirs(cache_dir, exist_ok=True)

    if gene_name.lower().strip() == "all":
        gene_key = "all"
        use_file = False
    else:
        genes = gene_name.split()
        genes = sorted(set(genes))
        gene_key_str = "\n".join(genes)  # one gene per line for hashing
        gene_key = hashlib.md5(gene_key_str.encode("utf-8")).hexdigest()
        use_file = True

    param_str = f"{reference}_{interval_expand}_{add_chr}"
    final_key_str = f"{gene_key}_{param_str}"
    final_hash = hashlib.md5(final_key_str.encode("utf-8")).hexdigest()

    cached_file = os.path.join(cache_dir, f"genes_{final_hash}.bed")

    if os.path.exists(cached_file):
        log_message("DEBUG", f"Found cached BED file: {cached_file}")
        return cached_file

    # If multiple genes, write them to a file
    gene_file_path = None
    if use_file:
        gene_file_path = os.path.join(cache_dir, f"genes_{final_hash}.txt")
        with open(gene_file_path, "w", encoding="utf-8") as gf:
            for g in genes:
                gf.write(g + "\n")

    # Generate BED
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "genes2bed", reference]
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])
    if use_file:
        # Use the file with genes
        cmd.extend(["-f", gene_file_path])

    log_message("INFO", f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, stdout=open(bed_path, "w", encoding="utf-8"), check=True)
    log_message("DEBUG", f"snpEff genes2bed completed, BED file at {bed_path}")

    # Sort the BED file
    sorted_bed = bed_path + ".sorted"
    sort_cmd = ["sortBed", "-i", bed_path]
    log_message("DEBUG", f"Sorting BED file with: {' '.join(sort_cmd)}")
    subprocess.run(sort_cmd, stdout=open(sorted_bed, "w", encoding="utf-8"), check=True)
    log_message("DEBUG", f"BED sorting completed, sorted file at {sorted_bed}")

    final_bed = sorted_bed
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
        final_bed = chr_bed
        log_message("DEBUG", f"Final BED file with 'chr': {final_bed}")

    shutil.move(final_bed, cached_file)
    log_message("DEBUG", f"Cached BED file created: {cached_file}")

    return cached_file

# File: variantcentrifuge/gene_bed.py
# Location: variantcentrifuge/variantcentrifuge/gene_bed.py

"""
Gene BED extraction module.

This module provides functionality to run snpEff genes2bed and produce
a BED file representing genes of interest. It also handles sorting and
optionally adding a 'chr' prefix to chromosome names.

Added caching logic:
- If genes != "all", sort them lexicographically, compute an MD5 hash,
  and use that to name the cached BED file.
- If "all", use a fixed hash like "all".
- Check if the cached file exists. If so, return it immediately.
- Otherwise, generate a new BED file, store it in the cache, and return it.
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
    If multiple genes are provided, they are sorted to create a stable cache key.

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

    # Prepare a cache directory for BED files
    cache_dir = os.path.join(output_dir, "bed_cache")
    os.makedirs(cache_dir, exist_ok=True)

    # Determine a stable hash for the gene list
    if gene_name.lower().strip() == "all":
        gene_key = "all"
    else:
        # Split by commas or whitespace, sort them
        genes = gene_name.replace(",", " ").split()
        genes = [g.strip() for g in genes if g.strip()]
        genes.sort()
        gene_key_str = " ".join(genes)
        # Compute MD5 hash of the sorted genes string
        gene_key = hashlib.md5(gene_key_str.encode("utf-8")).hexdigest()

    param_str = f"{reference}_{interval_expand}_{add_chr}"
    final_key_str = f"{gene_key}_{param_str}"
    final_hash = hashlib.md5(final_key_str.encode("utf-8")).hexdigest()

    # Cached file name
    cached_file = os.path.join(cache_dir, f"genes_{final_hash}.bed")

    # If cached file exists, just return it
    if os.path.exists(cached_file):
        log_message("DEBUG", f"Found cached BED file for these genes and params: {cached_file}")
        return cached_file

    # Otherwise, we need to generate the BED file
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "genes2bed", reference]
    if gene_name and gene_name.lower() != "all":
        cmd.append(gene_name)
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])

    log_message("INFO", f"Running: {' '.join(cmd)}")
    log_message("DEBUG", "Executing snpEff genes2bed command...")
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

    # Move the final_bed to the cached_file location using shutil.move
    shutil.move(final_bed, cached_file)
    log_message("DEBUG", f"Cached BED file created: {cached_file}")

    return cached_file

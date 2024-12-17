# variantcentrifuge/gene_bed.py

# This script normalizes gene inputs (single gene, multiple genes, or from a file)
# and uses snpEff genes2bed to generate a BED file representing the specified genes.
# It also supports caching to avoid regenerating BED files unnecessarily.

"""
Gene BED extraction and gene normalization module.

This module provides:
- normalize_genes: For normalizing gene inputs (single gene, multiple genes, or file).
- get_gene_bed: For generating a BED file corresponding to specified genes via snpEff genes2bed.
"""

import subprocess
import tempfile
import os
import hashlib
import shutil
import logging
import sys

logger = logging.getLogger("variantcentrifuge")


def normalize_genes(gene_name_str, gene_file_str, logger):
    """
    Normalize genes from either a single gene name, a list of genes,
    or a file containing gene names.

    If 'all' is provided or no genes after filtering, returns "all".

    Parameters
    ----------
    gene_name_str : str or None
        The gene name(s) provided via CLI.
    gene_file_str : str or None
        Path to a file containing gene names.
    logger : logging.Logger
        Logger for logging messages.

    Returns
    -------
    str
        A normalized, space-separated string of gene names or "all".
    """
    if gene_file_str and gene_file_str.strip():
        if not os.path.exists(gene_file_str):
            logger.error(f"Gene file {gene_file_str} not found.")
            sys.exit(1)
        genes_from_file = []
        with open(gene_file_str, "r", encoding="utf-8") as gf:
            for line in gf:
                line = line.strip()
                if line:
                    genes_from_file.append(line)
        genes = genes_from_file
    else:
        if not gene_name_str:
            logger.error("No gene name provided and no gene file provided.")
            sys.exit(1)
        # Prevent confusion if a file was incorrectly given to -g
        if os.path.exists(gene_name_str):
            logger.error(
                f"It looks like you provided a file '{gene_name_str}' to "
                f"-g/--gene-name."
            )
            logger.error(
                "If you meant to provide a file of gene names, please use "
                "-G/--gene-file instead."
            )
            sys.exit(1)

        g_str = gene_name_str.replace(",", " ")
        genes = [
            g.strip()
            for g_str_part in g_str.split()
            for g in [g_str_part.strip()] if g
        ]

    if len(genes) == 1 and genes[0].lower() == "all":
        return "all"

    # Sort and deduplicate
    genes = sorted(set(genes))
    if not genes:
        return "all"
    return " ".join(genes)


def get_gene_bed(reference, gene_name, interval_expand=0, add_chr=True, output_dir="output"):
    """
    Generate a BED file for the given gene(s) using snpEff genes2bed.
    If gene_name == "all", run once with no gene arguments.
    If multiple genes are provided, they are passed all at once as arguments.

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
    logger.debug(f"Entering get_gene_bed with reference={reference}, "
                 f"gene_name={gene_name}, interval_expand={interval_expand}, "
                 f"add_chr={add_chr}")

    cache_dir = os.path.join(output_dir, "bed_cache")
    os.makedirs(cache_dir, exist_ok=True)

    if gene_name.lower().strip() == "all":
        gene_key = "all"
        gene_args = []
    else:
        genes = gene_name.split()
        genes = sorted(genes)
        gene_key_str = " ".join(genes)
        gene_key = hashlib.md5(gene_key_str.encode("utf-8")).hexdigest()
        gene_args = genes

    param_str = f"{reference}_{interval_expand}_{add_chr}"
    final_key_str = f"{gene_key}_{param_str}"
    final_hash = hashlib.md5(final_key_str.encode("utf-8")).hexdigest()

    cached_file = os.path.join(cache_dir, f"genes_{final_hash}.bed")

    if os.path.exists(cached_file):
        logger.debug(f"Found cached BED file: {cached_file}")
        return cached_file

    # Generate BED
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "genes2bed", reference] + gene_args
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])

    logger.info(f"Running: {' '.join(cmd)}")
    with open(bed_path, "w", encoding="utf-8") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True, text=True)
    logger.debug(f"snpEff genes2bed completed, BED file at {bed_path}")

    # Sort the BED file
    sorted_bed = bed_path + ".sorted"
    sort_cmd = ["sortBed", "-i", bed_path]
    logger.debug(f"Sorting BED file with: {' '.join(sort_cmd)}")
    with open(sorted_bed, "w", encoding="utf-8") as out_f:
        subprocess.run(sort_cmd, stdout=out_f, check=True, text=True)
    logger.debug(f"BED sorting completed, sorted file at {sorted_bed}")

    final_bed = sorted_bed
    if add_chr:
        chr_bed = sorted_bed + ".chr"
        logger.debug(f"Adding 'chr' prefix to BED file {sorted_bed}")
        with open(chr_bed, "w", encoding="utf-8") as out_f, \
             open(sorted_bed, "r", encoding="utf-8") as in_f:
            for line in in_f:
                if not line.startswith("chr"):
                    out_f.write("chr" + line)
                else:
                    out_f.write(line)
        os.remove(sorted_bed)
        final_bed = chr_bed
        logger.debug(f"Final BED file with 'chr': {final_bed}")

    shutil.move(final_bed, cached_file)
    logger.debug(f"Cached BED file created: {cached_file}")

    return cached_file

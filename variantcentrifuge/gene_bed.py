# File: variantcentrifuge/gene_bed.py
# Location: variantcentrifuge/variantcentrifuge/gene_bed.py

"""
Gene BED extraction and gene normalization module.

This module provides:
- normalize_genes: For normalizing gene inputs (single gene, multiple genes, or file).
- get_gene_bed: For generating a BED file corresponding to specified genes via snpEff genes2bed.
"""

import hashlib
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from typing import List, Optional

logger = logging.getLogger("variantcentrifuge")


def normalize_genes(
    gene_name_str: Optional[str], gene_file_str: Optional[str], logger: logging.Logger
) -> str:
    """
    Normalize genes from either a single gene name, a list of genes,
    or a file containing gene names.

    If 'all' is provided or no genes after filtering, returns "all".

    Parameters
    ----------
    gene_name_str : str or None
        The gene name(s) provided via CLI (can be a single gene or space/comma-separated).
    gene_file_str : str or None
        Path to a file containing gene names, one per line.
    logger : logging.Logger
        Logger instance for logging messages.

    Returns
    -------
    str
        A normalized, space-separated string of gene names, or "all".

    Raises
    ------
    SystemExit
        If no gene name or file is provided, or if the specified file does not exist.
    """
    if gene_file_str and gene_file_str.strip():
        if not os.path.exists(gene_file_str):
            logger.error(f"Gene file {gene_file_str} not found.")
            sys.exit(1)
        genes_from_file: List[str] = []
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
        if os.path.exists(gene_name_str):
            logger.error(
                f"It looks like you provided a file '{gene_name_str}' to " f"-g/--gene-name."
            )
            logger.error(
                "If you meant to provide a file of gene names, please use "
                "-G/--gene-file instead."
            )
            sys.exit(1)

        g_str = gene_name_str.replace(",", " ")
        genes = [g.strip() for g_str_part in g_str.split() for g in [g_str_part.strip()] if g]

    if len(genes) == 1 and genes[0].lower() == "all":
        return "all"

    # Sort and deduplicate genes
    genes = sorted(set(genes))
    if not genes:
        return "all"
    return " ".join(genes)


def get_gene_bed(
    reference: str,
    gene_name: str,
    interval_expand: int = 0,
    add_chr: bool = True,
    output_dir: str = "output",
) -> str:
    """
    Generate a BED file for the given gene(s) using snpEff genes2bed.
    If gene_name == "all", the command runs without specifying genes.
    If multiple genes are provided, they are passed as arguments.

    Parameters
    ----------
    reference : str
        The reference genome name compatible with snpEff.
    gene_name : str
        "all" or space-separated list of gene names.
    interval_expand : int, optional
        Number of bases to expand upstream/downstream of the gene regions.
    add_chr : bool, optional
        Whether to add a 'chr' prefix to chromosome names in the BED file.
    output_dir : str, optional
        Directory to store cached BED files. Default is "output".

    Returns
    -------
    str
        Path to the final BED file.

    Raises
    ------
    subprocess.CalledProcessError
        If the snpEff genes2bed or sorting command fails.
    """
    logger.debug(
        "Entering get_gene_bed with reference=%s, gene_name=%s, interval_expand=%d, add_chr=%s",
        reference,
        gene_name,
        interval_expand,
        add_chr,
    )

    cache_dir = os.path.join(output_dir, "bed_cache")
    os.makedirs(cache_dir, exist_ok=True)

    if gene_name.lower().strip() == "all":
        gene_key = "all"
        gene_args: List[str] = []
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

    # Generate BED using snpEff genes2bed
    bed_fd, bed_path = tempfile.mkstemp(suffix=".bed")
    os.close(bed_fd)

    cmd = ["snpEff", "-Xmx8g", "genes2bed", reference] + gene_args
    if interval_expand > 0:
        cmd.extend(["-ud", str(interval_expand)])

    logger.info("Running: %s", " ".join(cmd))
    with open(bed_path, "w", encoding="utf-8") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True, text=True)
    logger.debug("snpEff genes2bed completed, BED file at %s", bed_path)

    # Sort the BED file
    sorted_bed = bed_path + ".sorted"
    sort_cmd = ["sortBed", "-i", bed_path]
    logger.debug("Sorting BED file with: %s", " ".join(sort_cmd))
    with open(sorted_bed, "w", encoding="utf-8") as out_f:
        subprocess.run(sort_cmd, stdout=out_f, check=True, text=True)
    logger.debug("BED sorting completed, sorted file at %s", sorted_bed)

    final_bed = sorted_bed
    if add_chr:
        chr_bed = sorted_bed + ".chr"
        logger.debug("Adding 'chr' prefix to BED file %s", sorted_bed)
        with open(chr_bed, "w", encoding="utf-8") as out_f, open(
            sorted_bed, "r", encoding="utf-8"
        ) as in_f:
            for line in in_f:
                if not line.startswith("chr"):
                    out_f.write("chr" + line)
                else:
                    out_f.write(line)
        os.remove(sorted_bed)
        final_bed = chr_bed
        logger.debug("Final BED file with 'chr': %s", final_bed)

    shutil.move(final_bed, cached_file)
    logger.debug("Cached BED file created: %s", cached_file)

    return cached_file

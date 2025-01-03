# File: variantcentrifuge/filters.py
# Location: variantcentrifuge/variantcentrifuge/filters.py

"""
Filtering module.

This module defines functions to extract variants from a VCF file using bcftools
and apply filters via SnpSift. Each function returns a filename containing the output.

Added:
- A function to filter final TSV rows by genotype (het, hom, comp_het), optionally
  applying different rules per gene based on a separate mapping file.
- Enhanced to append a reason label (e.g., "(comphet)") to each sample genotype if it
  passed the filter because of 'comp_het' or 'het' or 'hom'.

New changes for preserving extra fields (e.g., DP, AD) in parentheses:
- If a sample's entry is "325879(0/1:53,55:108)", we parse "0/1" as the "main" genotype
  but leave the rest of the substring (":53,55:108") intact. If this sample passes
  'het', we produce "325879(0/1:53,55:108)(het)" in the final output.
"""

import logging
import tempfile
import os
from typing import Dict, Any, Optional, Set, List
from .utils import run_command

logger = logging.getLogger("variantcentrifuge")


def extract_variants(
    vcf_file: str,
    bed_file: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Extract variants from a VCF using bcftools and a BED file, writing output
    to the specified compressed VCF ('.vcf.gz'). bcftools is invoked with the
    '-W' option, which writes the index file automatically.

    Parameters
    ----------
    vcf_file : str
        Path to the input VCF file.
    bed_file : str
        Path to the BED file containing genomic regions of interest.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.
        Expected keys include:
            - "threads": Number of threads to use with bcftools (default = 1).
    output_file : str
        Path to the final compressed output VCF file ('.vcf.gz').

    Returns
    -------
    str
        Path to the compressed VCF file (.vcf.gz) containing extracted variants.

    Raises
    ------
    RuntimeError
        If the extraction command fails.
    """
    threads = str(cfg.get("threads", 1))

    cmd = [
        "bcftools", "view",
        "--threads", threads,
        "-W",  # writes the index file automatically
        "-R", bed_file,
        "-Oz",  # compressed output
        "-o", output_file,
        vcf_file,
    ]
    logger.debug("Extracting variants with command: %s", " ".join(cmd))
    logger.debug("Output will be written to: %s", output_file)

    run_command(cmd)
    return output_file


def apply_snpsift_filter(
    variant_file: str,
    filter_string: str,
    cfg: Dict[str, Any],
    output_file: str
) -> str:
    """
    Apply a SnpSift filter to a variant file, then compress and index the output.

    Because our run_command function does not support shell pipelines,
    we split it into two steps:
      1) Write SnpSift filter output to a temporary uncompressed file (.vcf).
      2) Compress it with bgzip -@ <threads> to produce the final .vcf.gz.
      3) Index the resulting .vcf.gz with bcftools index.

    Parameters
    ----------
    variant_file : str
        Path to the compressed VCF file with extracted variants.
    filter_string : str
        SnpSift filter expression to apply.
    cfg : dict
        Configuration dictionary that may include paths and parameters for tools.
        Expected keys include:
            - "threads": Number of threads to use with bgzip and bcftools index (default = 1).
    output_file : str
        Path to the compressed VCF file (.vcf.gz) containing filtered variants.

    Returns
    -------
    str
        Path to the compressed VCF file (.vcf.gz) containing filtered variants.

    Raises
    ------
    RuntimeError
        If the filter command fails.
    """
    threads = str(cfg.get("threads", 1))

    # 1) Run SnpSift filter -> write uncompressed .vcf to a temporary file
    tmp_vcf = tempfile.mktemp(suffix=".vcf")
    snpsift_cmd = ["SnpSift", "filter", filter_string, variant_file]
    logger.debug("Applying SnpSift filter to produce uncompressed VCF: %s", tmp_vcf)
    run_command(snpsift_cmd, output_file=tmp_vcf)

    # 2) bgzip compress the temporary VCF
    bgzip_cmd = ["bgzip", "-@", threads, "-c", tmp_vcf]
    logger.debug("bgzip compressing to: %s", output_file)
    run_command(bgzip_cmd, output_file=output_file)

    # 3) bcftools index the resulting .vcf.gz
    index_cmd = ["bcftools", "index", "--threads", threads, output_file]
    logger.debug("Indexing filtered output with command: %s", " ".join(index_cmd))
    run_command(index_cmd)

    # Remove the uncompressed temp file
    if os.path.exists(tmp_vcf):
        os.remove(tmp_vcf)

    logger.debug("Filtered output is now available at: %s", output_file)
    return output_file


def filter_final_tsv_by_genotype(
    input_tsv: str,
    output_tsv: str,
    global_genotypes: Optional[Set[str]] = None,
    gene_genotype_file: Optional[str] = None,
    gene_column_name: str = "GENE",
    gt_column_name: str = "GT"
) -> None:
    """
    Filter the final TSV rows by genotype. This can be done globally (using a single
    set of requested genotypes like {"het"}, {"hom"}, or {"comp_het"}) or on a per-gene basis
    if gene_genotype_file is provided. The gene_genotype_file must contain at least two columns:
      1) GENE
      2) GENOTYPES (one or more of het, hom, comp_het, comma-separated)

    The logic is:
      - 'het'       => keep samples with genotype 0/1 or 1/0
      - 'hom'       => keep samples with genotype 1/1
      - 'comp_het'  => keep samples that have at least two distinct variants in the same gene
                       with genotype 0/1 or 1/0

    If a gene is defined in the gene_genotype_file, then the union of those genotype rules is applied.
    If a gene is not defined in that file (or if none is provided), global_genotypes is used.

    The resulting TSV keeps only lines that have at least one sample fulfilling the chosen
    genotype filters. If no samples remain on a line, that line is discarded.

    Additionally, if a sample passes because of 'het' or 'hom' or 'comp_het',
    we append a reason marker. For example:
      325879(0/1:53,55:108) => 325879(0/1:53,55:108)(het) or 325879(0/1:53,55:108)(het,comphet)

    Parameters
    ----------
    input_tsv : str
        Path to the input TSV, e.g. the final genotype_replaced TSV.
    output_tsv : str
        Path to the filtered output TSV.
    global_genotypes : set of str, optional
        A set of genotype filters to apply globally if no gene-specific rule is found.
        E.g. {"het"}, {"hom"}, or {"comp_het"}, or any combination thereof.
    gene_genotype_file : str, optional
        Path to a file with columns 'GENE' and 'GENOTYPES'. Each row can specify one or more
        genotype filters for a given gene, e.g. "BRCA2    het,comp_het".
    gene_column_name : str
        The column name in the TSV that specifies the gene.
    gt_column_name : str
        The column name in the TSV that specifies the genotype(s) for each sample.

    Returns
    -------
    None
        A new TSV is written to output_tsv.

    Notes
    -----
    For 'comp_het', we group rows by gene, identify samples that appear at least twice
    (with 0/1 or 1/0), then keep those rows for those samples only. We also annotate
    each genotype with the reason(s) it passed (het, hom, comphet).

    Example
    -------
    filter_final_tsv_by_genotype(
        "input.genotype_replaced.tsv",
        "output.genotype_filtered.tsv",
        global_genotypes={"comp_het"}
    )
    """
    if global_genotypes is None:
        global_genotypes = set()

    logger.debug("filter_final_tsv_by_genotype called with:")
    logger.debug("  input_tsv=%s", input_tsv)
    logger.debug("  output_tsv=%s", output_tsv)
    logger.debug("  global_genotypes=%s", global_genotypes)
    logger.debug("  gene_genotype_file=%s", gene_genotype_file)

    # Show the first few lines of the input_tsv (if available)
    try:
        with open(input_tsv, "r", encoding="utf-8") as inp_debug:
            logger.debug("First lines from input_tsv:")
            for i in range(3):
                line = inp_debug.readline()
                if not line:
                    break
                logger.debug("Line %d: %s", i + 1, line.rstrip("\n"))
    except FileNotFoundError:
        logger.error("input_tsv file not found: %s", input_tsv)
        raise

    # Read the gene -> genotype(s) mapping if provided
    gene_to_genotypes: Dict[str, Set[str]] = {}
    if gene_genotype_file and os.path.exists(gene_genotype_file):
        logger.debug("Attempting to read gene -> genotype rules from: %s", gene_genotype_file)
        with open(gene_genotype_file, "r", encoding="utf-8") as gfile:
            all_lines = gfile.readlines()

        logger.debug("First lines from gene_genotype_file:")
        for i, l in enumerate(all_lines[:3]):
            logger.debug("Line %d: %s", i + 1, l.rstrip("\n"))

        if not all_lines:
            logger.debug("Gene genotype file is empty, skipping parsing.")
        else:
            from io import StringIO
            gfile_replay = StringIO("".join(all_lines))
            header = next(gfile_replay).strip().split("\t")

            # Expect at least columns: GENE, GENOTYPES
            gene_idx = None
            geno_idx = None
            for i, col in enumerate(header):
                if col.upper() == "GENE":
                    gene_idx = i
                elif col.upper() == "GENOTYPES":
                    geno_idx = i
            if gene_idx is None or geno_idx is None:
                raise ValueError(
                    "gene_genotype_file must have columns named 'GENE' and 'GENOTYPES'."
                )
            line_count = 0
            for line in gfile_replay:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) <= max(gene_idx, geno_idx):
                    continue
                gname = parts[gene_idx].strip()
                genos = parts[geno_idx].replace(" ", "").split(",")
                if gname:
                    if gname not in gene_to_genotypes:
                        gene_to_genotypes[gname] = set()
                    gene_to_genotypes[gname].update(genos)
                line_count += 1
            logger.debug("Finished reading %d data lines from gene_genotype_file", line_count)
    else:
        logger.debug("No valid gene_genotype_file found, or file does not exist at: %s", gene_genotype_file)

    # Helpers to detect genotype as 'het' (0/1 or 1/0) or 'hom' (1/1)
    def is_het(gt_string: str) -> bool:
        return gt_string in ("0/1", "1/0")

    def is_hom(gt_string: str) -> bool:
        return gt_string == "1/1"

    # We'll parse the lines grouped by gene, so we can do 'comp_het' logic.
    lines_by_gene = {}

    with open(input_tsv, "r", encoding="utf-8") as inp:
        header = next(inp).rstrip("\n")
        header_cols = header.split("\t")
        # Identify gene and GT columns
        try:
            gene_idx = header_cols.index(gene_column_name)
        except ValueError:
            raise ValueError(
                f"Could not find gene column '{gene_column_name}' in TSV header."
            )
        try:
            gt_idx = header_cols.index(gt_column_name)
        except ValueError:
            raise ValueError(
                f"Could not find genotype column '{gt_column_name}' in TSV header."
            )

        for line in inp:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) <= max(gene_idx, gt_idx):
                continue

            gene_val = parts[gene_idx].strip()
            gt_val = parts[gt_idx].strip()

            # Parse each sample( genotype ) entry in GT
            sample_genotypes = {}
            if gt_val:
                entries = gt_val.split(";")
                for e in entries:
                    e = e.strip()
                    if not e:
                        continue
                    # example: "325879(0/1:53,55:108)"
                    # We'll store the entire substring "0/1:53,55:108" 
                    # but also parse out the "main genotype" before a colon if present
                    if "(" in e and ")" in e:
                        sample_name = e.split("(")[0].strip()
                        inside_paren = e.split("(")[1].strip(")")
                        sample_genotypes[sample_name] = inside_paren

            if gene_val not in lines_by_gene:
                lines_by_gene[gene_val] = []
            lines_by_gene[gene_val].append((line, sample_genotypes))

    # Identify which genes require 'comp_het'
    def gene_filters(g: str) -> Set[str]:
        if g in gene_to_genotypes:
            return gene_to_genotypes[g]
        return global_genotypes

    comp_het_qualified: Dict[str, Set[str]] = {}
    all_genes = set(lines_by_gene.keys())
    relevant_for_comp_het = {g for g in all_genes if "comp_het" in gene_filters(g)}

    # For each gene that uses comp_het, gather samples that have >=2 het variants
    for g in relevant_for_comp_het:
        sample_het_count = {}
        for line_str, sample_gt_dict in lines_by_gene[g]:
            for sample_name, genotype_substring in sample_gt_dict.items():
                # e.g. genotype_substring = "0/1:53,55:108"
                # parse the main genotype by splitting on the first colon
                main_gt = genotype_substring.split(":")[0] if ":" in genotype_substring else genotype_substring
                if is_het(main_gt):
                    sample_het_count[sample_name] = sample_het_count.get(sample_name, 0) + 1
        qualified_samples = {s for s, c in sample_het_count.items() if c >= 2}
        comp_het_qualified[g] = qualified_samples
        logger.debug("Gene '%s' => comp_het_qualified samples: %s", g, qualified_samples)

    # We'll build our filtered lines
    filtered_lines: List[str] = [header]  # keep original header as is

    # Now go through each gene's lines
    for g, items in lines_by_gene.items():
        filters_for_gene = gene_filters(g)
        do_het = "het" in filters_for_gene
        do_hom = "hom" in filters_for_gene
        do_comp_het = "comp_het" in filters_for_gene

        for line_str, sample_gt_dict in items:
            parts = line_str.split("\t")
            new_sample_entries = []

            # For each sample, check if it meets the filter(s)
            for sample_name, genotype_substring in sample_gt_dict.items():
                main_gt = genotype_substring.split(":")[0] if ":" in genotype_substring else genotype_substring

                reasons = []
                if do_het and is_het(main_gt):
                    reasons.append("het")
                if do_hom and is_hom(main_gt):
                    reasons.append("hom")
                if do_comp_het and is_het(main_gt) and sample_name in comp_het_qualified.get(g, set()):
                    reasons.append("comphet")

                if reasons:
                    # This sample passes the filter => we append the reason(s)
                    reason_str = f"({','.join(reasons)})"
                    # e.g. "325879(0/1:53,55:108)(het,comphet)"
                    new_sample_entries.append(f"{sample_name}({genotype_substring}){reason_str}")

            # If we have at least one sample that passes => keep the line
            if new_sample_entries:
                parts[gt_idx] = "; ".join(new_sample_entries)
                filtered_lines.append("\t".join(parts))

    # Write out the filtered lines
    with open(output_tsv, "w", encoding="utf-8") as out:
        for line in filtered_lines:
            out.write(line + "\n")

    logger.info(
        "Genotype filtering complete. Input: %s, Output: %s",
        input_tsv, output_tsv
    )

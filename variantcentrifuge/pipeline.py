# File: variantcentrifuge/pipeline.py
# Location: variantcentrifuge/variantcentrifuge/pipeline.py

"""
Pipeline orchestration module for variantcentrifuge.

This module provides high-level orchestration of the analysis steps:
- Checks external tools
- Normalizes genes
- Loads phenotypes and phenotype terms
- Extracts variants and applies filters
- Performs optional genotype replacement
- Integrates phenotype data
- Runs variant-level and gene-burden analyses using analyze_variants
- Generates metadata and optionally converts results to Excel
- Cleans up intermediates if requested

All steps are coordinated within the run_pipeline function.
"""

import argparse
import datetime
import gzip
import hashlib
import io
import logging
import locale
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union, Iterator

import pandas as pd

from .analyze_variants import analyze_variants
from .annotator import (
    load_custom_features,
    annotate_dataframe_with_features,
    validate_annotation_config,
)
from .converter import (
    append_tsv_as_sheet,
    convert_to_excel,
    finalize_excel_file,
    produce_report_json,
)
from .extractor import extract_fields
from .filters import (
    apply_snpsift_filter,
    extract_variants,
    filter_final_tsv_by_genotype,
    filter_tsv_with_expression,
)
from .gene_bed import get_gene_bed, normalize_genes
from .links import add_links_to_table
from .phenotype import aggregate_phenotypes_for_samples, load_phenotypes
from .phenotype_filter import filter_phenotypes
from .ped_reader import read_pedigree
from .replacer import replace_genotypes
from .scoring import read_scoring_config
from .utils import (
    check_external_tools,
    ensure_fields_in_extract,
    get_tool_version,
    normalize_snpeff_headers,
    run_command,
    sanitize_metadata_field,
    split_bed_file,
)
from variantcentrifuge.helpers import (
    annotate_variants_with_gene_lists,
    check_file,
    determine_case_control_sets,
    dump_df_to_xlsx,
    extract_gencode_id,
    get_vcf_names,
    get_vcf_regions,
    get_vcf_samples,
    get_vcf_size,
    match_IGV_link_columns,
    read_sequencing_manifest,
)

# Import the SNPeff annotation splitting function
from .vcf_eff_one_per_line import process_vcf_file

logger = logging.getLogger("variantcentrifuge")


def smart_open(filename: str, mode: str = "r", encoding: str = "utf-8"):
    """
    Open a file with automatic gzip support based on file extension.

    Parameters
    ----------
    filename : str
        Path to the file
    mode : str
        File opening mode ('r', 'w', 'rt', 'wt', etc.)
    encoding : str
        Text encoding (for text modes)

    Returns
    -------
    file object
        Opened file handle
    """
    if filename.endswith(".gz"):
        # Ensure text mode for gzip
        if "t" not in mode and "b" not in mode:
            mode = mode + "t"
        return gzip.open(filename, mode, encoding=encoding)
    else:
        # For regular files, only add encoding for text mode
        if "b" not in mode:
            return open(filename, mode, encoding=encoding)
        else:
            return open(filename, mode)


def sort_tsv_by_gene(
    input_file: str,
    output_file: str,
    gene_column: str = "GENE",
    temp_dir: Optional[str] = None,
    memory_limit: str = "2G",
    parallel: int = 4,
) -> str:
    """
    Sort a TSV file by gene column using external sort command for memory efficiency.

    This function uses the system sort command which is memory-efficient and can
    handle files larger than available RAM by using disk-based sorting.

    Parameters
    ----------
    input_file : str
        Path to the input TSV file (can be gzipped)
    output_file : str
        Path to the output sorted TSV file (can be gzipped)
    gene_column : str
        Name of the gene column to sort by
    temp_dir : str, optional
        Directory for temporary files during sorting
    memory_limit : str
        Memory limit for sort command (e.g., '2G', '500M')
    parallel : int
        Number of parallel threads for sorting

    Returns
    -------
    str
        Path to the sorted output file
    """
    logger.info(f"Sorting TSV file by gene column: {input_file} -> {output_file}")
    logger.info(f"Using memory limit: {memory_limit}, parallel threads: {parallel}")

    # First, we need to find the gene column index
    with smart_open(input_file, "r") as f:
        header = f.readline().strip()
        if not header:
            logger.error(f"Input file '{input_file}' is empty or has no header")
            raise ValueError(f"Input file '{input_file}' is empty or has no header")
        columns = header.split("\t")

    try:
        gene_col_idx = columns.index(gene_column) + 1  # sort uses 1-based indexing
    except ValueError:
        logger.error(f"Gene column '{gene_column}' not found in TSV header")
        raise ValueError(f"Gene column '{gene_column}' not found in TSV file")

    logger.info(f"Found gene column '{gene_column}' at position {gene_col_idx}")

    # Build sort arguments with memory efficiency options
    sort_args = [
        f"-k{gene_col_idx},{gene_col_idx}",  # Sort by gene column
        f"--buffer-size={memory_limit}",  # Memory limit
        f"--parallel={parallel}",  # Parallel threads
        "--stable",  # Stable sort for consistent results
        "--compress-program=gzip",  # Use gzip for temp files
    ]

    if temp_dir:
        sort_args.append(f"-T {temp_dir}")

    sort_cmd = " ".join(sort_args)

    # Escape file paths to prevent shell injection
    import shlex

    safe_input = shlex.quote(input_file)
    safe_output = shlex.quote(output_file)

    # Build the complete command
    if input_file.endswith(".gz"):
        if output_file.endswith(".gz"):
            # Both gzipped
            cmd = (
                f"gzip -cd {safe_input} | "
                f"{{ IFS= read -r header; echo \"$header\"; sort {sort_cmd} -t$'\\t'; }} | "
                f"gzip -c > {safe_output}"
            )
        else:
            # Input gzipped, output not
            cmd = (
                f"gzip -cd {safe_input} | "
                f"{{ IFS= read -r header; echo \"$header\"; sort {sort_cmd} -t$'\\t'; }} "
                f"> {safe_output}"
            )
    else:
        if output_file.endswith(".gz"):
            # Input not gzipped, output gzipped
            cmd = (
                f'{{ IFS= read -r header < {safe_input}; echo "$header"; tail -n +2 {safe_input} | '
                f"sort {sort_cmd} -t$'\\t'; }} | "
                f"gzip -c > {safe_output}"
            )
        else:
            # Neither gzipped
            cmd = (
                f'{{ IFS= read -r header < {safe_input}; echo "$header"; tail -n +2 {safe_input} | '
                f"sort {sort_cmd} -t$'\\t'; }} "
                f"> {safe_output}"
            )

    # Execute the command
    logger.debug(f"Running sort command: {cmd}")

    import subprocess
    import shutil

    # Find bash executable
    bash_path = shutil.which("bash") or "/bin/bash"
    if not os.path.exists(bash_path):
        logger.warning(f"bash not found at {bash_path}, falling back to shell default")
        bash_path = None

    # Use bash explicitly to ensure bash-specific syntax works
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, executable=bash_path)

    if result.returncode != 0:
        logger.error(f"Sort command failed: {result.stderr}")
        raise RuntimeError(f"Failed to sort TSV file: {result.stderr}")

    # Verify output file was created
    if not os.path.exists(output_file):
        logger.error(f"Sort command completed but output file '{output_file}' was not created")
        raise RuntimeError(f"Output file '{output_file}' was not created")

    # Check if output file has content
    with smart_open(output_file, "r") as f:
        first_line = f.readline()
        if not first_line:
            logger.error(f"Sort command completed but output file '{output_file}' is empty")
            raise RuntimeError(f"Output file '{output_file}' is empty")

    logger.info("Successfully sorted TSV file by gene column")
    return output_file


def process_chunked_pipeline(
    final_tsv: str,
    final_output: str,
    cfg: Dict[str, Any],
    custom_features: Dict[str, Any],
    scoring_config: Optional[Dict[str, Any]],
    pedigree_data: Optional[Dict[str, Any]],
    args: argparse.Namespace,
    base_name: str,
    intermediate_dir: str,
    chunksize: int = 10000,
) -> None:
    """
    Process the pipeline in a streaming, chunked manner to handle large files.

    This function replaces the in-memory DataFrame processing with a streaming
    approach that processes data in gene-aware chunks.

    Parameters
    ----------
    final_tsv : str
        Path to the input TSV file
    final_output : str
        Path to the final output file
    cfg : dict
        Configuration dictionary
    custom_features : dict
        Custom annotation features
    scoring_config : dict, optional
        Scoring configuration
    pedigree_data : dict, optional
        Pedigree data for inheritance analysis
    args : argparse.Namespace
        Command line arguments
    base_name : str
        Base name for output files
    intermediate_dir : str
        Directory for intermediate files
    chunksize : int
        Target chunk size for processing
    """
    logger.info("Starting chunked pipeline processing")

    # Determine gene column name
    with smart_open(final_tsv, "r") as f:
        header_line = f.readline().strip()
        header_cols = header_line.split("\t")

    gene_col_name = None
    for col in header_cols:
        if col == "GENE" or col.endswith(".GENE"):
            gene_col_name = col
            break

    if not gene_col_name:
        logger.warning("No GENE column found, falling back to non-chunked processing")
        # Fall back to original processing
        return False

    # Prepare output file
    output_handle = smart_open(final_output, "w")
    header_written = False

    # Prepare analyze_variants config
    temp_cfg = cfg.copy()
    temp_cfg["perform_gene_burden"] = False
    temp_cfg["scoring_config"] = scoring_config
    temp_cfg["pedigree_data"] = pedigree_data
    temp_cfg["calculate_inheritance"] = cfg.get("calculate_inheritance", False)
    temp_cfg["sample_list"] = cfg.get("sample_list", "")

    # Process chunks
    total_chunks = 0
    total_variants = 0

    try:
        for chunk_num, chunk_df in enumerate(
            read_tsv_in_gene_chunks(final_tsv, gene_column=gene_col_name, chunksize=chunksize)
        ):
            logger.info(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")

            # Apply unified custom annotations
            if custom_features and (
                bool(custom_features.get("regions_by_chrom"))
                or bool(custom_features.get("gene_lists"))
                or bool(custom_features.get("json_gene_data"))
            ):
                chunk_df = annotate_dataframe_with_features(chunk_df, custom_features)
            else:
                if "Custom_Annotation" not in chunk_df.columns:
                    chunk_df["Custom_Annotation"] = ""

            # Convert chunk to temporary file for analyze_variants
            chunk_tsv = os.path.join(intermediate_dir, f"{base_name}.chunk_{chunk_num}.tsv.gz")
            chunk_df.to_csv(chunk_tsv, sep="\t", index=False, na_rep="", compression="gzip")

            # Run analyze_variants on chunk
            chunk_buffer = []
            with smart_open(chunk_tsv, "r") as inp:
                for line in analyze_variants(inp, temp_cfg):
                    chunk_buffer.append(line)

            # Clean up chunk file
            if not cfg.get("keep_intermediates"):
                os.remove(chunk_tsv)

            # Process chunk results
            if chunk_buffer:
                # Apply late filtering if needed
                if (
                    cfg.get("late_filtering", False)
                    and cfg.get("filters")
                    and len(chunk_buffer) > 1
                ):
                    # Convert to DataFrame for filtering
                    chunk_str = "\n".join(chunk_buffer)
                    filter_df = pd.read_csv(
                        io.StringIO(chunk_str), sep="\t", dtype=str, keep_default_na=False
                    )

                    # Apply filter using the existing function
                    from .filters import filter_dataframe_with_query

                    if hasattr(
                        filter_dataframe_with_query, "__name__"
                    ):  # Check if it's the pandas query version
                        filter_df = filter_dataframe_with_query(filter_df, cfg["filters"])

                    # Convert back to buffer
                    output_str = io.StringIO()
                    filter_df.to_csv(output_str, sep="\t", index=False, na_rep="")
                    chunk_buffer = output_str.getvalue().rstrip("\n").split("\n")

                # Add links if needed
                if not cfg.get("no_links", False) and len(chunk_buffer) > 1:
                    links_config = cfg.get("links", {})
                    if links_config:
                        chunk_buffer = add_links_to_table(chunk_buffer, links_config)

                # Process inheritance output if needed
                if (
                    cfg.get("calculate_inheritance", False)
                    and cfg.get("inheritance_mode")
                    and len(chunk_buffer) > 1
                ):
                    chunk_str = "\n".join(chunk_buffer)
                    inheritance_df = pd.read_csv(
                        io.StringIO(chunk_str), sep="\t", dtype=str, keep_default_na=False
                    )

                    from .inheritance.analyzer import process_inheritance_output

                    inheritance_df = process_inheritance_output(
                        inheritance_df, cfg.get("inheritance_mode", "simple")
                    )

                    output_str = io.StringIO()
                    inheritance_df.to_csv(output_str, sep="\t", index=False, na_rep="")
                    chunk_buffer = output_str.getvalue().rstrip("\n").split("\n")

                # Write chunk results
                if not header_written and chunk_buffer:
                    # First chunk - process header for VAR_ID
                    header_line = chunk_buffer[0]
                    header_parts = header_line.split("\t")

                    # Add any additional columns requested
                    if hasattr(args, "add_column") and args.add_column:
                        header_parts.extend(args.add_column)

                    # Add VAR_ID as first column
                    new_header = ["VAR_ID"] + header_parts
                    output_handle.write("\t".join(new_header) + "\n")
                    header_written = True

                    # Find column indices for VAR_ID generation
                    header_indices = {col: idx for idx, col in enumerate(header_parts)}
                    chrom_idx = header_indices.get("CHROM", None)
                    pos_idx = header_indices.get("POS", None)
                    ref_idx = header_indices.get("REF", None)
                    alt_idx = header_indices.get("ALT", None)

                    # Write data rows with VAR_ID
                    for line in chunk_buffer[1:]:
                        if line.strip():
                            fields = line.split("\t")

                            # Generate VAR_ID
                            chrom_val = fields[chrom_idx] if chrom_idx is not None else ""
                            pos_val = fields[pos_idx] if pos_idx is not None else ""
                            ref_val = fields[ref_idx] if ref_idx is not None else ""
                            alt_val = fields[alt_idx] if alt_idx is not None else ""

                            combined = f"{chrom_val}{pos_val}{ref_val}{alt_val}"
                            short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]
                            var_id = f"var_{total_variants + 1:04d}_{short_hash}"

                            # Add any additional blank columns if requested
                            if hasattr(args, "add_column") and args.add_column:
                                fields.extend([""] * len(args.add_column))

                            # Write line with VAR_ID
                            new_fields = [var_id] + fields
                            output_handle.write("\t".join(new_fields) + "\n")
                            total_variants += 1
                else:
                    # Subsequent chunks - skip header but add VAR_ID
                    for line in chunk_buffer[1:]:
                        if line.strip():
                            fields = line.split("\t")

                            # Generate VAR_ID (using the already established indices)
                            chrom_val = fields[chrom_idx] if chrom_idx is not None else ""
                            pos_val = fields[pos_idx] if pos_idx is not None else ""
                            ref_val = fields[ref_idx] if ref_idx is not None else ""
                            alt_val = fields[alt_idx] if alt_idx is not None else ""

                            combined = f"{chrom_val}{pos_val}{ref_val}{alt_val}"
                            short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]
                            var_id = f"var_{total_variants + 1:04d}_{short_hash}"

                            # Add any additional blank columns if requested
                            if hasattr(args, "add_column") and args.add_column:
                                fields.extend([""] * len(args.add_column))

                            # Write line with VAR_ID
                            new_fields = [var_id] + fields
                            output_handle.write("\t".join(new_fields) + "\n")
                            total_variants += 1

            total_chunks += 1

            # Log progress
            if total_chunks % 10 == 0:
                logger.info(f"Processed {total_chunks} chunks, {total_variants} variants written")

    finally:
        output_handle.close()

    logger.info(
        f"Chunked processing complete: {total_chunks} chunks, {total_variants} total variants"
    )

    # Apply final filter if provided
    if cfg.get("final_filter"):
        logger.info("Applying final filter to results")

        # Read the output file
        df = pd.read_csv(
            final_output,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            compression="gzip" if final_output.endswith(".gz") else None,
        )

        # Apply the filter
        from .filters import filter_dataframe_with_query

        df = filter_dataframe_with_query(df, cfg["final_filter"])

        # Write back
        df.to_csv(
            final_output,
            sep="\t",
            index=False,
            na_rep="",
            compression="gzip" if final_output.endswith(".gz") else None,
        )

        logger.info(f"Final filter retained {len(df)} variants")

    return True


def read_tsv_in_gene_chunks(
    filepath: str, gene_column: str = "GENE", chunksize: int = 10000, compression: str = "infer"
) -> Iterator[pd.DataFrame]:
    """
    Read a TSV file in gene-aware chunks, ensuring all variants for a gene are in the same chunk.

    This generator function reads a TSV file chunk by chunk but ensures that all rows
    belonging to the same gene are yielded together. This is critical for analyses
    like compound heterozygous detection that require all variants for a gene to be
    processed together.

    Parameters
    ----------
    filepath : str
        Path to the TSV file (can be gzipped)
    gene_column : str
        Name of the column containing gene identifiers
    chunksize : int
        Target number of rows per chunk (actual chunks may be larger to keep genes together)
    compression : str
        Compression type ('infer' to auto-detect, 'gzip', or None)

    Yields
    ------
    pd.DataFrame
        DataFrames containing complete gene data, never splitting a gene across chunks

    Notes
    -----
    The input TSV file MUST be sorted by the gene column for this to work correctly.
    Each yielded chunk will contain at least one complete gene and may contain multiple
    genes if they fit within the chunksize limit.
    """
    logger.info(f"Starting gene-aware chunked reading of {filepath}")
    logger.info(f"Target chunk size: {chunksize} rows, gene column: {gene_column}")

    # Initialize variables
    gene_buffer = pd.DataFrame()
    chunks_yielded = 0
    total_rows = 0

    # Use pandas chunked reader
    chunk_iterator = pd.read_csv(
        filepath,
        sep="\t",
        chunksize=chunksize,
        dtype=str,
        keep_default_na=False,
        compression=compression,
    )

    for chunk_num, chunk in enumerate(chunk_iterator):
        # Add chunk to buffer
        if gene_buffer.empty:
            gene_buffer = chunk
        else:
            gene_buffer = pd.concat([gene_buffer, chunk], ignore_index=True)

        # Warn if buffer is getting very large (potential memory issue)
        if len(gene_buffer) > chunksize * 10:
            logger.warning(
                f"Gene buffer has grown to {len(gene_buffer)} rows. "
                f"Consider increasing chunksize if you have genes with many variants."
            )

        # Check if we have enough data to yield
        while len(gene_buffer) >= chunksize:
            # Find the last complete gene in the buffer
            if gene_column not in gene_buffer.columns:
                logger.error(f"Gene column '{gene_column}' not found in TSV file")
                raise ValueError(f"Gene column '{gene_column}' not found in TSV file")

            # Get unique genes in order of appearance
            gene_values = gene_buffer[gene_column].values

            # Find where gene changes occur
            gene_change_indices = [0]  # Start with first index
            current_gene = gene_values[0]

            for i in range(1, len(gene_values)):
                if gene_values[i] != current_gene:
                    gene_change_indices.append(i)
                    current_gene = gene_values[i]

            # Find the split point - we want at least chunksize rows but complete genes
            split_index = None
            for idx in gene_change_indices:
                if idx >= chunksize:
                    split_index = idx
                    break

            # If no split point found (all rows belong to same gene or last gene extends beyond chunksize)
            if split_index is None:
                # Keep accumulating - we'll yield this in the final cleanup
                break
            else:
                # Yield the complete genes
                yield_df = gene_buffer.iloc[:split_index].copy()
                yield yield_df

                chunks_yielded += 1
                total_rows += len(yield_df)

                # Keep the remaining data
                gene_buffer = gene_buffer.iloc[split_index:].reset_index(drop=True)

                # Log progress
                if chunks_yielded % 10 == 0:
                    logger.info(f"Yielded {chunks_yielded} chunks, processed {total_rows} rows")

    # Yield any remaining data
    if not gene_buffer.empty:
        yield gene_buffer
        chunks_yielded += 1
        total_rows += len(gene_buffer)

    logger.info(
        f"Completed gene-aware chunked reading: {chunks_yielded} chunks, {total_rows} total rows"
    )


def remove_vcf_extensions(filename: str) -> str:
    """
    Remove common VCF-related extensions from a filename.

    Parameters
    ----------
    filename : str
        The input filename, possibly ending in .vcf, .vcf.gz, or .gz.

    Returns
    -------
    str
        The filename base without VCF-related extensions.
    """
    if filename.endswith(".vcf.gz"):
        return filename[:-7]
    elif filename.endswith(".vcf"):
        return filename[:-4]
    elif filename.endswith(".gz"):
        return filename[:-3]
    return filename


def compute_base_name(vcf_path: str, gene_name: str) -> str:
    """
    Compute a base name for output files based on the VCF filename and genes.

    If multiple genes are specified, create a hash to represent them.
    If 'all' is specified, append '.all'.
    Otherwise, append the gene name if it's not already in the VCF base name.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.
    gene_name : str
        The normalized gene name string.

    Returns
    -------
    str
        A base name for output files.
    """
    genes = gene_name.strip()
    vcf_base = os.path.basename(vcf_path)
    vcf_base = remove_vcf_extensions(vcf_base)

    if genes.lower() == "all":
        return f"{vcf_base}.all"
    split_genes = genes.split()
    if len(split_genes) > 1:
        gene_hash = hashlib.md5(genes.encode("utf-8")).hexdigest()[:8]
        return f"{vcf_base}.multiple-genes-{gene_hash}"
    else:
        if split_genes and split_genes[0].lower() in vcf_base.lower():
            return vcf_base
        else:
            return f"{vcf_base}.{split_genes[0]}" if split_genes else vcf_base


def load_terms_from_file(file_path: Optional[str], logger: logging.Logger) -> List[str]:
    """
    Load terms (HPO terms, sample IDs, etc.) from a file, one per line.

    Parameters
    ----------
    file_path : str or None
        Path to a file containing one term per line.
    logger : logging.Logger
        Logger instance for error logging.

    Returns
    -------
    list of str
        A list of terms loaded from the file.

    Raises
    ------
    SystemExit
        If the file is missing or empty and a file_path was specified.
    """
    terms: List[str] = []
    if file_path:
        if not os.path.exists(file_path):
            logger.error(f"Required file not found: {file_path}")
            sys.exit(1)
        with open(file_path, "r", encoding="utf-8") as f:
            found_any = False
            for line in f:
                t = line.strip()
                if t:
                    found_any = True
                    terms.append(t)
            if not found_any:
                logger.error(f"File {file_path} is empty or invalid.")
                sys.exit(1)
    return terms


def parse_samples_from_vcf(vcf_file: str) -> List[str]:
    """
    Parse sample names from a VCF file by reading its header line.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.

    Returns
    -------
    list of str
        A list of sample names extracted from the VCF header.
    """
    output = run_command(["bcftools", "view", "-h", vcf_file], output_file=None)
    chrom_line = None
    for line in output.splitlines():
        if line.startswith("#CHROM"):
            chrom_line = line
            break
    if not chrom_line:
        return []
    fields = chrom_line.strip().split("\t")
    if len(fields) <= 9:
        return []
    samples = fields[9:]
    return samples


def _process_bed_chunk(
    chunk_index: int,
    chunk_bed_file: str,
    original_vcf: str,
    base_name: str,
    intermediate_dir: str,
    cfg: Dict[str, Any],
    transcripts: List[str] = None,
    args: argparse.Namespace = None,
) -> str:
    """Worker function to process a single BED chunk."""
    chunk_base_name = f"{base_name}.chunk_{chunk_index}"
    worker_log_prefix = f"[Worker {chunk_index}]"
    logger.info(
        f"{worker_log_prefix} Starting processing for BED chunk: {os.path.basename(chunk_bed_file)}"
    )

    # Define chunk-specific intermediate file paths
    chunk_variants_vcf = os.path.join(intermediate_dir, f"{chunk_base_name}.variants.vcf.gz")
    chunk_splitted_before_vcf = os.path.join(
        intermediate_dir, f"{chunk_base_name}.splitted_before_filters.vcf.gz"
    )
    chunk_splitted_after_vcf = os.path.join(
        intermediate_dir, f"{chunk_base_name}.splitted_after_filters.vcf.gz"
    )
    chunk_filtered_vcf = os.path.join(intermediate_dir, f"{chunk_base_name}.filtered.vcf.gz")
    chunk_transcript_filtered_vcf = os.path.join(
        intermediate_dir, f"{chunk_base_name}.transcript_filtered.vcf.gz"
    )
    chunk_extracted_tsv = os.path.join(intermediate_dir, f"{chunk_base_name}.extracted.tsv.gz")

    # Run the standard filtering and extraction pipeline on this chunk
    # Use limited threads for bcftools since it only helps with compression
    chunk_cfg = cfg.copy()
    chunk_cfg["threads"] = min(2, cfg.get("threads", 1))  # Max 2 threads for bcftools
    # Extract variants (with optional bcftools pre-filter)
    extract_variants(original_vcf, chunk_bed_file, chunk_cfg, chunk_variants_vcf)

    # Handle snpeff splitting mode
    splitting_mode = cfg.get("snpeff_splitting_mode", None)

    if splitting_mode == "before_filters":
        logger.info(
            f"{worker_log_prefix} Splitting multiple SNPeff (EFF/ANN) annotations before main filtering."
        )
        from .vcf_eff_one_per_line import process_vcf_file

        process_vcf_file(chunk_variants_vcf, chunk_splitted_before_vcf)
        prefiltered_for_snpsift = chunk_splitted_before_vcf
    else:
        prefiltered_for_snpsift = chunk_variants_vcf

    # Apply filter only if late filtering is not enabled
    if cfg.get("late_filtering", False):
        logger.info(
            f"{worker_log_prefix} Late filtering enabled - skipping early SnpSift filter step"
        )
        import shutil

        shutil.copy2(prefiltered_for_snpsift, chunk_filtered_vcf)
    else:
        apply_snpsift_filter(prefiltered_for_snpsift, cfg["filters"], chunk_cfg, chunk_filtered_vcf)

    # If splitting after filters
    if splitting_mode == "after_filters":
        logger.info(
            f"{worker_log_prefix} Splitting multiple SNPeff (EFF/ANN) annotations after main filters."
        )
        from .vcf_eff_one_per_line import process_vcf_file

        process_vcf_file(chunk_filtered_vcf, chunk_splitted_after_vcf)
        post_filter_for_transcripts = chunk_splitted_after_vcf
    else:
        post_filter_for_transcripts = chunk_filtered_vcf

    # Apply transcript filter if transcripts are provided
    if transcripts:
        logger.info(f"{worker_log_prefix} Filtering for transcripts using SnpSift.")
        or_clauses = [f"(EFF[*].TRID = '{tr}')" for tr in transcripts]
        transcript_filter_expr = " | ".join(or_clauses)

        apply_snpsift_filter(
            post_filter_for_transcripts,
            transcript_filter_expr,
            chunk_cfg,
            chunk_transcript_filtered_vcf,
        )
        final_filtered_for_extraction = chunk_transcript_filtered_vcf
    else:
        final_filtered_for_extraction = post_filter_for_transcripts

    # Extract fields
    # If user wants to append extra sample fields, ensure they're in the main fields
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields", []):
        updated_fields = ensure_fields_in_extract(
            cfg["fields_to_extract"], cfg["extra_sample_fields"]
        )
        field_list = " ".join(updated_fields.strip().split())
    else:
        field_list = " ".join(
            (cfg["fields_to_extract"] or (args.fields if args else "")).strip().split()
        )

    logger.debug(f"{worker_log_prefix} Extracting fields: {field_list}")

    snpsift_sep = cfg.get("extract_fields_separator", ":")
    chunk_cfg["extract_fields_separator"] = snpsift_sep
    extract_fields(final_filtered_for_extraction, field_list, chunk_cfg, chunk_extracted_tsv)

    logger.info(f"{worker_log_prefix} Finished. Output TSV at {chunk_extracted_tsv}")
    return chunk_extracted_tsv


def _run_parallel_pipeline(
    args, cfg, gene_name, num_workers, base_name, intermediate_dir, transcripts=None
):
    """Handle the parallel execution of the pipeline by splitting the master BED file."""
    logger.info("Generating master BED file for all specified genes...")
    master_bed_file = get_gene_bed(
        cfg["reference"],
        gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True),
        output_dir=args.output_dir,
    )

    if not os.path.exists(master_bed_file) or os.path.getsize(master_bed_file) == 0:
        logger.error("Gene BED file could not be created or is empty. Check genes or reference.")
        sys.exit(1)

    # Check if requested genes found in reference
    if gene_name.lower() != "all":
        requested_genes = gene_name.split()
        found_genes = set()
        with open(master_bed_file, "r", encoding="utf-8") as bf:
            for line in bf:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    g_name = parts[3].split(";")[0].strip()
                    found_genes.add(g_name)
        missing = [g for g in requested_genes if g not in found_genes]
        if missing:
            logger.warning(
                "The following gene(s) were not found in the reference: " f"{', '.join(missing)}"
            )
            if cfg.get("debug_level", "INFO") == "ERROR":
                sys.exit(1)

    logger.info(f"Splitting master BED file into {num_workers} chunks...")
    chunk_bed_files = split_bed_file(master_bed_file, num_workers, intermediate_dir)

    if not chunk_bed_files:
        logger.error("Failed to split BED file into chunks. Falling back to single-threaded mode.")
        return _run_single_thread_pipeline(
            args, cfg, gene_name, base_name, intermediate_dir, transcripts
        )

    chunk_tsv_files = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {
            executor.submit(
                _process_bed_chunk,
                i,
                bed_chunk,
                args.vcf_file,
                base_name,
                intermediate_dir,
                cfg,
                transcripts,
                args,
            ): i
            for i, bed_chunk in enumerate(chunk_bed_files)
        }
        for future in as_completed(futures):
            try:
                result_path = future.result()
                chunk_tsv_files.append(result_path)
            except Exception as e:
                logger.error(f"A worker process failed: {e}")
                executor.shutdown(wait=False, cancel_futures=True)
                sys.exit("Pipeline aborted due to worker failure.")

    combined_tsv_path = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv.gz")
    logger.info(
        f"All chunks processed. Merging {len(chunk_tsv_files)} files into {combined_tsv_path}..."
    )

    # Merge chunks into combined file
    try:
        with smart_open(combined_tsv_path, "w") as outfile:
            first_file = True
            for chunk_file in sorted(chunk_tsv_files):
                # Verify chunk file exists before trying to read
                if not os.path.exists(chunk_file):
                    logger.error(f"Chunk file {chunk_file} not found during merge")
                    raise FileNotFoundError(f"Chunk file {chunk_file} not found")

                with smart_open(chunk_file, "r") as infile:
                    if first_file:
                        import shutil

                        shutil.copyfileobj(infile, outfile)
                        first_file = False
                    else:
                        next(infile)  # Skip header
                        shutil.copyfileobj(infile, outfile)
    except Exception as e:
        logger.error(f"Failed to merge chunk files: {e}")
        raise

    # Cleanup intermediate files only after successful merge
    if not cfg.get("keep_intermediates"):
        logger.debug("Cleaning up intermediate chunk files")
        for file_path in chunk_tsv_files + chunk_bed_files:
            try:
                if os.path.exists(file_path):
                    os.remove(file_path)
                    logger.debug(f"Removed intermediate file: {file_path}")
            except OSError as e:
                logger.warning(f"Could not remove intermediate file {file_path}: {e}")

    return combined_tsv_path


def _run_single_thread_pipeline(
    args, cfg, gene_name, base_name, intermediate_dir, transcripts=None
):
    """Run the original linear pipeline logic."""
    # This function is now incorporated into the main run_pipeline function
    # when use_parallel is False. This stub is kept for compatibility.
    raise NotImplementedError("This function should not be called directly anymore.")


def archive_results_directory(output_dir: str, base_name: str) -> Optional[str]:
    """
    Archive the results directory into a compressed tar.gz file.
    
    Parameters
    ----------
    output_dir : str
        Path to the directory to archive
    base_name : str
        Base name for the archive (from input VCF)
        
    Returns
    -------
    Optional[str]
        Path to the created archive file, or None if archiving failed
    """
    import tarfile
    from datetime import datetime
    
    # Generate timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create archive name: variantcentrifuge_results_<basename>_<timestamp>.tar.gz
    archive_name = f"variantcentrifuge_results_{base_name}_{timestamp}.tar.gz"
    
    # Place archive in parent directory of output_dir
    parent_dir = os.path.dirname(os.path.abspath(output_dir))
    archive_path = os.path.join(parent_dir, archive_name)
    
    try:
        logger.info(f"Creating archive: {archive_path}")
        
        with tarfile.open(archive_path, "w:gz") as tar:
            # Add the entire output directory
            tar.add(output_dir, arcname=os.path.basename(output_dir))
            
        # Verify archive was created and get size
        if os.path.exists(archive_path):
            size_mb = os.path.getsize(archive_path) / (1024 * 1024)
            logger.info(f"Archive created successfully: {archive_path} ({size_mb:.1f} MB)")
            return archive_path
        else:
            logger.error("Archive file was not created")
            return None
            
    except Exception as e:
        logger.error(f"Failed to create archive: {str(e)}")
        return None


def run_pipeline(
    args: argparse.Namespace, cfg: Dict[str, Any], start_time: datetime.datetime
) -> None:
    """
    High-level orchestration of the pipeline steps.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    cfg : dict
        Configuration dictionary merged from CLI and config file.
    start_time : datetime.datetime
        The start time of the run.

    Returns
    -------
    None
        Writes output files and may print results to stdout.
    """
    # Check external tools
    check_external_tools()

    # Normalize genes
    gene_name = normalize_genes(args.gene_name if args.gene_name else "", args.gene_file, logger)
    logger.debug(f"Normalized gene list: {gene_name}")

    # Determine base name and output paths
    if args.output_file is not None:
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = os.path.join(args.output_dir, os.path.basename(args.output_file))
    else:
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, gene_name)
        final_output = os.path.join(args.output_dir, f"{base_name}.final.tsv")

    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    if not cfg["no_stats"] and not args.stats_output_file:
        cfg["stats_output_file"] = os.path.join(intermediate_dir, f"{base_name}.statistics.tsv")
    else:
        cfg["stats_output_file"] = args.stats_output_file

    # Load phenotypes if provided
    phenotypes = {}
    use_phenotypes = False
    if args.phenotype_file and args.phenotype_sample_column and args.phenotype_value_column:
        phenotypes = load_phenotypes(
            args.phenotype_file,
            args.phenotype_sample_column,
            args.phenotype_value_column,
        )
        if not phenotypes:
            logger.error(
                f"No phenotype data loaded from {args.phenotype_file}. "
                "Check file formatting and column names."
            )
            sys.exit(1)
        use_phenotypes = True

    # Store phenotypes in cfg for determine_case_control_sets usage
    cfg["phenotypes"] = phenotypes

    # Load scoring configuration if path is provided
    scoring_config = None
    if cfg.get("scoring_config_path"):
        try:
            scoring_config = read_scoring_config(cfg["scoring_config_path"])
            logger.info(
                f"Successfully loaded scoring configuration from {cfg['scoring_config_path']}"
            )
        except Exception as e:
            logger.error(f"Failed to load scoring configuration: {e}")
            sys.exit(1)

    # Load pedigree data if provided
    pedigree_data = None
    if cfg.get("calculate_inheritance"):
        if cfg.get("ped_file"):
            try:
                pedigree_data = read_pedigree(cfg["ped_file"])
                cfg["pedigree_data"] = pedigree_data
                logger.info(f"Loaded pedigree data for {len(pedigree_data)} individuals.")
            except Exception as e:
                logger.error(f"Failed to load PED file: {e}")
                sys.exit(1)
        else:
            # Create single-sample pedigree data for inheritance analysis
            logger.info("No PED file provided - treating all samples as affected individuals.")
            # When no PED file is provided, create empty pedigree_data
            # Will be populated later when we have sample information
            pedigree_data = {}
            cfg["pedigree_data"] = pedigree_data

    # Validate and load custom annotation features
    annotation_errors = validate_annotation_config(cfg)
    if annotation_errors:
        logger.error("Custom annotation configuration errors:")
        for error in annotation_errors:
            logger.error(f"  - {error}")
        sys.exit(1)

    custom_features = load_custom_features(cfg)

    # Load phenotype terms
    case_hpo_terms = []
    control_hpo_terms = []
    if args.case_phenotypes:
        case_hpo_terms = [t.strip() for t in args.case_phenotypes.split(",") if t.strip()]
    case_hpo_terms += load_terms_from_file(args.case_phenotypes_file, logger)

    if args.control_phenotypes:
        control_hpo_terms = [t.strip() for t in args.control_phenotypes.split(",") if t.strip()]
    control_hpo_terms += load_terms_from_file(args.control_phenotypes_file, logger)

    cfg["case_phenotypes"] = case_hpo_terms
    cfg["control_phenotypes"] = control_hpo_terms

    # Load explicit sample sets
    case_samples = []
    control_samples = []
    if args.case_samples:
        case_samples = [s.strip() for s in args.case_samples.split(",") if s.strip()]
    case_samples += load_terms_from_file(args.case_samples_file, logger)

    if args.control_samples:
        control_samples = [s.strip() for s in args.control_samples.split(",") if s.strip()]
    control_samples += load_terms_from_file(args.control_samples_file, logger)

    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    # Determine whether to use parallel or single-threaded processing
    num_workers = cfg.get("threads", 1)
    use_parallel = num_workers > 1

    if use_parallel:
        logger.info(f"Running in parallel mode with {num_workers} workers.")
    else:
        logger.info("Running in single-threaded mode.")

    # Parse transcripts before running pipeline (needed for both parallel and single-threaded)
    transcripts = []
    if args.transcript_list:
        transcripts.extend([t.strip() for t in args.transcript_list.split(",") if t.strip()])
    if args.transcript_file:
        if not os.path.exists(args.transcript_file):
            logger.error(f"Transcript file not found: {args.transcript_file}")
            sys.exit(1)
        with open(args.transcript_file, "r", encoding="utf-8") as tf:
            for line in tf:
                tr = line.strip()
                if tr:
                    transcripts.append(tr)
    transcripts = list(set(transcripts))  # remove duplicates
    if transcripts:
        logger.info(f"Will filter for {len(transcripts)} transcript(s).")

    # Pre-configure fields to extract (needed for parallel workers)
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields", []):
        updated_fields = ensure_fields_in_extract(
            cfg["fields_to_extract"], cfg["extra_sample_fields"]
        )
        cfg["fields_to_extract"] = updated_fields
        logger.debug(
            f"Updated fields_to_extract with extra sample fields: {cfg['fields_to_extract']}"
        )

    # Filenames for intermediate steps
    variants_file = os.path.join(intermediate_dir, f"{base_name}.variants.vcf.gz")
    splitted_before_file = os.path.join(
        intermediate_dir, f"{base_name}.splitted_before_filters.vcf.gz"
    )
    splitted_after_file = os.path.join(
        intermediate_dir, f"{base_name}.splitted_after_filters.vcf.gz"
    )
    filtered_file = os.path.join(intermediate_dir, f"{base_name}.filtered.vcf.gz")
    transcript_filtered_file = os.path.join(
        intermediate_dir, f"{base_name}.transcript_filtered.vcf.gz"
    )
    extracted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv.gz")
    genotype_replaced_tsv = os.path.join(intermediate_dir, f"{base_name}.genotype_replaced.tsv.gz")
    phenotype_added_tsv = os.path.join(intermediate_dir, f"{base_name}.phenotypes_added.tsv.gz")
    gene_burden_tsv = os.path.join(args.output_dir, f"{base_name}.gene_burden.tsv")

    # Parse samples from VCF
    original_samples = parse_samples_from_vcf(args.vcf_file)
    if cfg.get("remove_sample_substring"):
        substring_to_remove = cfg["remove_sample_substring"]
        if substring_to_remove and substring_to_remove.strip():
            original_samples = [s.replace(substring_to_remove, "") for s in original_samples]

    # If no PED file was provided, populate pedigree_data with samples as affected
    if pedigree_data is not None and not pedigree_data and cfg.get("calculate_inheritance", False):
        logger.debug("Populating pedigree data for single sample analysis")
        for sample in original_samples:
            pedigree_data[sample] = {
                "sample_id": sample,
                "father_id": "0",
                "mother_id": "0",
                "sex": "0",  # Unknown sex
                "affected_status": "2",  # Affected
            }

    # -----------------------------------------------------------------------
    # Step 1: Run parallel or single-threaded pipeline to get extracted TSV
    # -----------------------------------------------------------------------
    if use_parallel:
        # Parallel pipeline processing
        extracted_tsv = _run_parallel_pipeline(
            args, cfg, gene_name, num_workers, base_name, intermediate_dir, transcripts
        )
    else:
        # Single-threaded pipeline processing (original behavior)
        # Extract gene BED
        bed_file = get_gene_bed(
            cfg["reference"],
            gene_name,
            interval_expand=cfg.get("interval_expand", 0),
            add_chr=cfg.get("add_chr", True),
            output_dir=args.output_dir,
        )
        logger.debug(f"Gene BED created at: {bed_file}")

        if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
            logger.error(
                "Gene BED file could not be created or is empty. Check genes or reference."
            )
            sys.exit(1)

        # Check if requested genes found in reference
        if gene_name.lower() != "all":
            requested_genes = gene_name.split()
            found_genes = set()
            with open(bed_file, "r", encoding="utf-8") as bf:
                for line in bf:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        g_name = parts[3].split(";")[0].strip()
                        found_genes.add(g_name)
            missing = [g for g in requested_genes if g not in found_genes]
            if missing:
                logger.warning(
                    "The following gene(s) were not found in the reference: "
                    f"{', '.join(missing)}"
                )
                if cfg.get("debug_level", "INFO") == "ERROR":
                    sys.exit(1)

        # Extract variants (with optional bcftools pre-filter)
        extract_variants(args.vcf_file, bed_file, cfg, variants_file)

        # Handle snpeff splitting mode (None, 'before_filters', or 'after_filters')
        splitting_mode = cfg.get("snpeff_splitting_mode", None)

        if splitting_mode == "before_filters":
            logger.info("Splitting multiple SNPeff (EFF/ANN) annotations before main filtering.")
            process_vcf_file(variants_file, splitted_before_file)
            prefiltered_for_snpsift = splitted_before_file
        else:
            # either None or 'after_filters'
            prefiltered_for_snpsift = variants_file

        # Apply SnpSift filter => filtered_file (skip if late filtering is enabled)
        if cfg.get("late_filtering", False):
            logger.info("Late filtering enabled - skipping early SnpSift filter step")
            # Just copy the file without filtering
            import shutil

            shutil.copy2(prefiltered_for_snpsift, filtered_file)
        else:
            apply_snpsift_filter(prefiltered_for_snpsift, cfg["filters"], cfg, filtered_file)

        # If splitting after filters
        if splitting_mode == "after_filters":
            logger.info(
                "Splitting multiple SNPeff (EFF/ANN) annotations after main filters, before transcript filter."
            )
            process_vcf_file(filtered_file, splitted_after_file)
            post_filter_for_transcripts = splitted_after_file
        else:
            post_filter_for_transcripts = filtered_file

        # Transcript filtering

        if transcripts:
            logger.info("Filtering for transcripts using SnpSift.")
            or_clauses = [f"(EFF[*].TRID = '{tr}')" for tr in transcripts]
            transcript_filter_expr = " | ".join(or_clauses)

            apply_snpsift_filter(
                post_filter_for_transcripts,
                transcript_filter_expr,
                cfg,
                transcript_filtered_file,
            )
            final_filtered_for_extraction = transcript_filtered_file
        else:
            final_filtered_for_extraction = post_filter_for_transcripts

        # Extract fields => extracted_tsv
        field_list = " ".join((cfg["fields_to_extract"] or args.fields).strip().split())
        logger.debug(f"Extracting fields: {field_list} -> {extracted_tsv}")

        snpsift_sep = cfg.get("extract_fields_separator", ":")
        logger.debug(f"Using SnpSift subfield separator for extractFields: '{snpsift_sep}'")

        temp_cfg = cfg.copy()
        temp_cfg["extract_fields_separator"] = snpsift_sep
        extract_fields(final_filtered_for_extraction, field_list, temp_cfg, extracted_tsv)

    # Get field_list for sorting logic (needed for both parallel and single-threaded)
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields", []):
        updated_fields = ensure_fields_in_extract(
            cfg["fields_to_extract"], cfg["extra_sample_fields"]
        )
        field_list_for_sorting = " ".join(updated_fields.strip().split())
    else:
        field_list_for_sorting = " ".join((cfg["fields_to_extract"] or args.fields).strip().split())

    # Sort the extracted TSV by gene for efficient chunked processing
    sorted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.sorted.tsv.gz")

    # Check if GENE column exists in the extracted fields (could be GENE or ANN[0].GENE etc)
    possible_gene_fields = ["GENE", "ANN[0].GENE", "ANN[*].GENE"]
    gene_field_found = None
    for field in possible_gene_fields:
        if field in field_list_for_sorting:
            gene_field_found = field
            break

    if gene_field_found:
        # Check actual column name in the file (SnpSift may normalize it)
        with smart_open(extracted_tsv, "r") as f:
            header_line = f.readline().strip()
            header_cols = header_line.split("\t")

        # Find the actual gene column name
        gene_col_name = None
        for col in header_cols:
            if col == "GENE" or col.endswith(".GENE"):
                gene_col_name = col
                break

        if gene_col_name:
            logger.info(
                f"Sorting extracted TSV by gene column '{gene_col_name}' for efficient processing..."
            )
            sort_tsv_by_gene(
                extracted_tsv,
                sorted_tsv,
                gene_column=gene_col_name,
                temp_dir=intermediate_dir,
                memory_limit=cfg.get("sort_memory_limit", "2G"),
                parallel=cfg.get("sort_parallel", 4),
            )
            # Replace extracted_tsv with sorted version
            if not cfg.get("keep_intermediates"):
                os.remove(extracted_tsv)
            extracted_tsv = sorted_tsv
        else:
            logger.warning("Could not find GENE column in TSV header, skipping sorting step")
    else:
        logger.warning(
            "GENE field not found in extracted fields configuration, skipping sorting step"
        )

    # Check if GT column is present
    with smart_open(extracted_tsv, "r") as f:
        header_line = f.readline().strip()
    columns = header_line.split("\t")
    gt_present = "GT" in columns

    cfg["sample_list"] = ",".join(original_samples) if original_samples else ""

    # MODIFIED: Perform inheritance analysis BEFORE genotype replacement
    if cfg.get("calculate_inheritance", False) and pedigree_data is not None and original_samples:
        logger.info("Performing inheritance analysis before genotype replacement...")

        try:
            # Load the extracted TSV with individual sample genotypes
            import pandas as pd

            df = pd.read_csv(
                extracted_tsv,
                sep="\t",
                dtype=str,
                keep_default_na=False,
                compression="gzip" if extracted_tsv.endswith(".gz") else None,
            )

            # Check if we have a single GT column with colon-separated values (from GEN[*].GT)
            # or individual sample columns
            if "GT" in df.columns and len(df) > 0:
                # Check if GT column contains colon-separated genotypes
                first_gt = str(df["GT"].iloc[0]) if not df["GT"].isna().all() else ""
                snpsift_sep = cfg.get("extract_fields_separator", ":")

                if snpsift_sep in first_gt:
                    # GT column contains colon-separated genotypes for all samples
                    # We need to split it into individual sample columns
                    logger.debug(
                        f"GT column contains {snpsift_sep}-separated values, splitting into sample columns"
                    )

                    # Pre-create all sample columns data
                    sample_data = {sample_id: [] for sample_id in original_samples}

                    # Extract genotypes for each row
                    for idx, row in df.iterrows():
                        gt_value = str(row["GT"])
                        if gt_value and gt_value != "NA" and gt_value != "nan":
                            genotypes = gt_value.split(snpsift_sep)
                            if len(genotypes) != len(original_samples):
                                logger.warning(
                                    f"Row {idx}: Expected {len(original_samples)} genotypes but found {len(genotypes)}"
                                )
                            for i, sample_id in enumerate(original_samples):
                                if i < len(genotypes):
                                    sample_data[sample_id].append(genotypes[i])
                                else:
                                    sample_data[sample_id].append("./.")
                        else:
                            # Missing GT data
                            for sample_id in original_samples:
                                sample_data[sample_id].append("./.")

                    # Create a new DataFrame with sample columns and concatenate
                    sample_df = pd.DataFrame(sample_data)
                    df = pd.concat([df, sample_df], axis=1)
                else:
                    # GT column doesn't contain separator, might be indexed columns
                    logger.debug(
                        "GT column doesn't contain separator, checking for indexed columns"
                    )
                    rename_map = {}
                    for i, sample_id in enumerate(original_samples):
                        if f"GT_{i}" in df.columns:
                            rename_map[f"GT_{i}"] = sample_id
                    if rename_map:
                        logger.debug(f"Renaming indexed columns: {list(rename_map.keys())[:5]}")
                        df = df.rename(columns=rename_map)
            else:
                # No GT column, check for indexed columns
                rename_map = {}
                for i, sample_id in enumerate(original_samples):
                    if f"GT_{i}" in df.columns:
                        rename_map[f"GT_{i}"] = sample_id
                if rename_map:
                    logger.debug(
                        f"Found indexed columns without GT column: {list(rename_map.keys())[:5]}"
                    )
                    df = df.rename(columns=rename_map)

            # Perform inheritance analysis
            from .inheritance.analyzer import analyze_inheritance

            use_vectorized = not args.no_vectorized_comp_het
            df = analyze_inheritance(
                df, pedigree_data, original_samples, use_vectorized_comp_het=use_vectorized
            )

            # Remove individual sample columns before saving
            # But keep all other columns including case/control counts
            cols_to_keep = [col for col in df.columns if col not in original_samples]
            df = df[cols_to_keep]

            # Save the TSV with inheritance columns added
            df.to_csv(
                extracted_tsv,
                sep="\t",
                index=False,
                na_rep="",
                compression="gzip" if extracted_tsv.endswith(".gz") else None,
            )
            logger.info("Inheritance analysis complete. Results added to extracted TSV.")

            # Set flag to skip inheritance analysis later
            cfg["inheritance_already_calculated"] = True

        except Exception as e:
            logger.error(f"Failed to perform inheritance analysis: {e}")
            logger.error("Inheritance analysis will be skipped")
            cfg["calculate_inheritance"] = False

    # Genotype replacement logic
    if not args.no_replacement and gt_present:
        logger.info("Starting genotype replacement...")
        replacement_start_time = time.time()
        lines_written = 0

        try:
            with smart_open(extracted_tsv, "r") as inp, smart_open(
                genotype_replaced_tsv, "w"
            ) as out:
                for line in replace_genotypes(inp, cfg):
                    out.write(line + "\n")
                    lines_written += 1

                    # Log progress every 10000 lines
                    if lines_written % 10000 == 0:
                        elapsed = time.time() - replacement_start_time
                        rate = lines_written / elapsed if elapsed > 0 else 0
                        logger.info(
                            f"Genotype replacement: {lines_written} lines in {elapsed:.1f}s ({rate:.0f} lines/sec)"
                        )
                        out.flush()  # Ensure data is written to disk
        except Exception as e:
            logger.error(f"Error during genotype replacement at line {lines_written}: {e}")
            logger.error("Consider using --no-replacement flag to skip this step")
            raise

        total_time = time.time() - replacement_start_time
        logger.info(f"Genotype replacement completed: {lines_written} lines in {total_time:.1f}s")
        replaced_tsv = genotype_replaced_tsv
    else:
        replaced_tsv = extracted_tsv

    # If user has appended extra fields, they might want them removed after genotype assembly
    if cfg.get("append_extra_sample_fields", False) and cfg.get("extra_sample_fields"):
        logger.debug(
            "User config => removing columns after replacement: %s",
            cfg["extra_sample_fields"],
        )
        stripped_tsv = os.path.join(intermediate_dir, f"{base_name}.stripped_extras.tsv.gz")

        with smart_open(replaced_tsv, "r") as inp, smart_open(stripped_tsv, "w") as out:
            # Read the header line from replaced_tsv
            raw_header_line = next(inp).rstrip("\n")
            original_header_cols = raw_header_line.split("\t")

            # Create a normalized copy of the header to see how columns might have been renamed
            # by SnpSift or other tools:
            # Keep normalized header code for potential future use
            # normed_header_line = normalize_snpeff_headers([raw_header_line])[0]

            # Build a map from normalized_col -> original_index
            normalized_to_index = {}
            for i, col in enumerate(original_header_cols):
                # Single-column approach to re-using normalize_snpeff_headers:
                normed_col = normalize_snpeff_headers([col])[0]
                normalized_to_index[normed_col] = i

            remove_indices = []
            # For each user-supplied field, normalize it and see if it matches a normalized header col
            for raw_col_name in cfg["extra_sample_fields"]:
                single_line = normalize_snpeff_headers([raw_col_name])[0]
                if single_line in normalized_to_index:
                    idx = normalized_to_index[single_line]
                    remove_indices.append(idx)
                    logger.debug(
                        "Removing normalized column '%s' => real index %d",
                        single_line,
                        idx,
                    )
                else:
                    logger.warning(
                        "Column '%s' was requested for removal but not found in header!",
                        raw_col_name,
                    )

            remove_indices.sort(reverse=True)
            new_header = [h for i, h in enumerate(original_header_cols) if i not in remove_indices]
            out.write("\t".join(new_header) + "\n")

            for line_num, line in enumerate(inp, start=2):
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line + "\n")
                    continue
                parts = line.split("\t")
                for idx in remove_indices:
                    if idx < len(parts):
                        parts.pop(idx)
                out.write("\t".join(parts) + "\n")

        os.rename(stripped_tsv, replaced_tsv)
    logger.debug("Extra column removal (if requested) complete.")

    # Integrate phenotypes if provided
    if use_phenotypes:
        pattern = re.compile(r"^([^()]+)(?:\([^)]+\))?$")

        with smart_open(replaced_tsv, "r") as inp, smart_open(phenotype_added_tsv, "w") as out:
            header = next(inp).rstrip("\n")
            header_fields = header.split("\t")
            header_fields.append("phenotypes")
            out.write("\t".join(header_fields) + "\n")

            wrote_data = False
            gt_idx = header_fields.index("GT") if "GT" in header_fields else None
            for line in inp:
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line + "\n")
                    continue
                fields_line = line.split("\t")

                if gt_idx is not None and gt_idx < len(fields_line):
                    gt_value = fields_line[gt_idx]
                    samples_in_line = []
                    if gt_value.strip():
                        sample_entries = gt_value.split(";")
                        for entry in sample_entries:
                            entry = entry.strip()
                            if entry:
                                m = pattern.match(entry)
                                if m:
                                    s = m.group(1).strip()
                                    if s:
                                        samples_in_line.append(s)
                    pheno_str = ""
                    if samples_in_line:
                        pheno_str = aggregate_phenotypes_for_samples(samples_in_line, phenotypes)
                    fields_line.append(pheno_str)
                else:
                    fields_line.append("")
                out.write("\t".join(fields_line) + "\n")
                wrote_data = True

            if not wrote_data:
                logger.warning(
                    "Phenotype integration found no data rows to process. Only header row will be written."
                )
        final_tsv = phenotype_added_tsv
    else:
        final_tsv = replaced_tsv

    # Genotype filtering if requested
    if getattr(args, "genotype_filter", None) or getattr(args, "gene_genotype_file", None):
        genotype_filtered_tsv = os.path.join(
            args.output_dir, f"{base_name}.genotype_filtered.tsv.gz"
        )
        genotype_modes = set()
        if getattr(args, "genotype_filter", None):
            genotype_modes = set(g.strip() for g in args.genotype_filter.split(",") if g.strip())
        filter_final_tsv_by_genotype(
            input_tsv=final_tsv,
            output_tsv=genotype_filtered_tsv,
            global_genotypes=genotype_modes,
            gene_genotype_file=args.gene_genotype_file,
        )
        final_tsv = genotype_filtered_tsv

    # If Excel requested
    if args.xlsx:
        excel_file = os.path.splitext(final_output or f"{base_name}.final.tsv")[0] + ".xlsx"
        cfg["final_excel_file"] = excel_file
    else:
        cfg["final_excel_file"] = None

    # Determine if we should use chunked processing
    use_chunked_processing = not cfg.get("no_chunked_processing", False)

    # Check file size to decide on processing mode
    if use_chunked_processing:
        file_size_mb = os.path.getsize(final_tsv) / (1024 * 1024)
        logger.info(f"Input file size: {file_size_mb:.1f} MB")

        # Use chunked processing for files > 100MB or if explicitly requested
        if file_size_mb > 100 or cfg.get("force_chunked_processing", False):
            logger.info("Using chunked processing for large file")

            # Prepare intermediate output file
            chunked_output = os.path.join(intermediate_dir, f"{base_name}.chunked_output.tsv.gz")

            # Run chunked processing
            success = process_chunked_pipeline(
                final_tsv=final_tsv,
                final_output=chunked_output,
                cfg=cfg,
                custom_features=custom_features,
                scoring_config=scoring_config,
                pedigree_data=pedigree_data,
                args=args,
                base_name=base_name,
                intermediate_dir=intermediate_dir,
                chunksize=cfg.get("chunk_size", 10000),
            )

            if success:
                # Read the chunked output into buffer for final processing
                with smart_open(chunked_output, "r") as f:
                    buffer = [line.rstrip("\n") for line in f]

                # Clean up intermediate file
                if not cfg.get("keep_intermediates"):
                    os.remove(chunked_output)
            else:
                # Fall back to regular processing
                use_chunked_processing = False
        else:
            logger.info("File size below threshold, using regular processing")
            use_chunked_processing = False

    # Regular (non-chunked) processing
    if not use_chunked_processing:
        # Load TSV into DataFrame for unified annotation and analysis
        import pandas as pd

        try:
            df = pd.read_csv(
                final_tsv,
                sep="\t",
                dtype=str,
                keep_default_na=False,
                compression="gzip" if final_tsv.endswith(".gz") else None,
            )
            logger.info(f"Loaded {len(df)} variants for annotation and analysis")
        except Exception as e:
            logger.error(f"Failed to load TSV file {final_tsv}: {e}")
            sys.exit(1)

        # Apply unified custom annotations
        has_custom_features = (
            bool(custom_features.get("regions_by_chrom"))
            or bool(custom_features.get("gene_lists"))
            or bool(custom_features.get("json_gene_data"))
        )

        if has_custom_features:
            logger.info("Applying unified custom annotations...")
            df = annotate_dataframe_with_features(df, custom_features)
        else:
            # Add empty Custom_Annotation column for consistency
            df["Custom_Annotation"] = ""

        # Inheritance analysis already done before genotype replacement if requested
        # No need to re-augment from VCF or re-calculate

        # Convert DataFrame back to file format for analyze_variants
        annotated_tsv = os.path.join(intermediate_dir, f"{base_name}.annotated.tsv.gz")
        df.to_csv(annotated_tsv, sep="\t", index=False, na_rep="", compression="gzip")

        # Run analyze_variants for variant-level results
        temp_cfg = cfg.copy()
        temp_cfg["perform_gene_burden"] = False
        temp_cfg["scoring_config"] = scoring_config  # Pass the loaded config
        temp_cfg["pedigree_data"] = pedigree_data  # Pass the loaded pedigree data
        temp_cfg["calculate_inheritance"] = cfg.get("calculate_inheritance", False)
        temp_cfg["sample_list"] = cfg.get("sample_list", "")  # Ensure sample_list is passed
        buffer = []
        with smart_open(annotated_tsv, "r") as inp:
            line_count = 0
            for line in analyze_variants(inp, temp_cfg):
                buffer.append(line)
                line_count += 1
            if line_count <= 1:  # Only header or nothing
                logger.warning(
                    "No variant-level results produced after filtering. Generating empty output files with headers."
                )

    # Note: Gene list annotation is now handled by the unified annotation system above
    # This ensures --annotate-gene-list files are processed through the new Custom_Annotation column
    # The old separate column system is deprecated but kept here for backward compatibility checking
    gene_list_annotation_files = cfg.get("annotate_gene_list_files", [])
    if gene_list_annotation_files and len(buffer) > 1:
        logger.warning(
            "Gene list annotation is now integrated into the unified Custom_Annotation system. "
            "Individual gene list columns are no longer added separately."
        )

    # Apply late filtering if enabled and filters are specified
    if cfg.get("late_filtering", False) and cfg.get("filters") and len(buffer) > 1:
        logger.info("Applying late filtering on scored and annotated data...")

        # Write buffer to a temporary file
        temp_scored_tsv = os.path.join(intermediate_dir, f"{base_name}.scored.tsv.gz")
        with smart_open(temp_scored_tsv, "w") as f:
            for line in buffer:
                f.write(line + "\n")

        # Apply the filter
        temp_filtered_tsv = os.path.join(intermediate_dir, f"{base_name}.late_filtered.tsv.gz")
        filter_tsv_with_expression(
            temp_scored_tsv,
            temp_filtered_tsv,
            cfg["filters"],
            pandas_query=False,  # Use SnpSift-style syntax
        )

        # Read the filtered results back into buffer
        with smart_open(temp_filtered_tsv, "r") as f:
            buffer = [line.rstrip("\n") for line in f]

        if len(buffer) <= 1:
            logger.warning("No variants passed late filtering")

    # Add links if not disabled and we have data rows
    if not cfg.get("no_links", False) and len(buffer) > 1:
        links_config = cfg.get("links", {})
        if links_config:
            logger.debug("Adding link columns to final output.")
            buffer = add_links_to_table(buffer, links_config)
    else:
        logger.debug("Link columns are disabled by configuration or no variants to link.")

    # Insert a leading VAR_ID column
    def add_variant_identifier(lines: List[str]) -> List[str]:
        """
        Insert a leading VAR_ID column into each line.

        The VAR_ID is generated from a zero-padded incrementing counter plus
        a short 4-character hash based on CHROM, POS, REF, and ALT.

        Parameters
        ----------
        lines : list of str
            The lines of the final variant table.

        Returns
        -------
        list of str
            A new list of lines with a VAR_ID column inserted at the beginning.
        """
        new_lines = []
        header_indices = {}
        row_counter = 1

        # Identify header and relevant columns
        header_line = lines[0].rstrip("\n")
        header_parts = header_line.split("\t")
        for idx, col in enumerate(header_parts):
            header_indices[col] = idx

        # Prepend header
        new_header = ["VAR_ID"] + header_parts
        new_lines.append("\t".join(new_header))

        chrom_idx = header_indices.get("CHROM", None)
        pos_idx = header_indices.get("POS", None)
        ref_idx = header_indices.get("REF", None)
        alt_idx = header_indices.get("ALT", None)

        for line in lines[1:]:
            line = line.rstrip("\n")
            if not line.strip():
                new_lines.append(line)
                continue

            fields = line.split("\t")
            chrom_val = fields[chrom_idx] if chrom_idx is not None else ""
            pos_val = fields[pos_idx] if pos_idx is not None else ""
            ref_val = fields[ref_idx] if ref_idx is not None else ""
            alt_val = fields[alt_idx] if alt_idx is not None else ""

            combined = f"{chrom_val}{pos_val}{ref_val}{alt_val}"
            short_hash = hashlib.md5(combined.encode("utf-8")).hexdigest()[:4]

            var_id = f"var_{row_counter:04d}_{short_hash}"
            row_counter += 1
            new_line = [var_id] + fields
            new_lines.append("\t".join(new_line))

        return new_lines

    # Skip post-processing if chunked processing was used (already handled)
    chunked_processing_used = False
    if use_chunked_processing and "buffer" in locals() and buffer:
        # Check if buffer already has VAR_ID column (indicates chunked processing completed)
        if buffer[0].startswith("VAR_ID"):
            chunked_processing_used = True

    if not chunked_processing_used:
        # Process inheritance output based on mode if inheritance was calculated
        if (
            cfg.get("calculate_inheritance", False)
            and cfg.get("inheritance_mode")
            and len(buffer) > 1
        ):
            # Convert buffer to DataFrame for processing
            import io

            buffer_str = "\n".join(buffer)
            df_for_inheritance = pd.read_csv(
                io.StringIO(buffer_str), sep="\t", dtype=str, keep_default_na=False
            )

            # Apply inheritance output processing
            from .inheritance.analyzer import process_inheritance_output

            df_for_inheritance = process_inheritance_output(
                df_for_inheritance, cfg.get("inheritance_mode", "simple")
            )

            # Convert back to buffer
            output_str = io.StringIO()
            df_for_inheritance.to_csv(output_str, sep="\t", index=False, na_rep="")
            buffer = output_str.getvalue().rstrip("\n").split("\n")

        # Add variant identifiers only if we have data rows
        if len(buffer) > 1:
            buffer = add_variant_identifier(buffer)
        else:
            # For empty results, just add the VAR_ID column to the header
            if buffer and "\t" in buffer[0]:
                header_parts = buffer[0].split("\t")
                new_header = ["VAR_ID"] + header_parts
                buffer[0] = "\t".join(new_header)

    # Append named blank columns if requested
    def add_named_columns(lines: List[str], col_names: List[str]) -> List[str]:
        """
        Append columns named in col_names to the table, with blank data for each row.

        Parameters
        ----------
        lines : list of str
            Lines (header + data) of the final table.
        col_names : list of str
            The column names to add.

        Returns
        -------
        list of str
            The updated table lines with extra named blank columns.
        """
        if not lines or not col_names:
            return lines

        new_lines = []
        header_line = lines[0].rstrip("\n")
        header_parts = header_line.split("\t")
        header_parts.extend(col_names)
        new_header = "\t".join(header_parts)
        new_lines.append(new_header)

        for line in lines[1:]:
            line = line.rstrip("\n")
            if not line.strip():
                new_lines.append(line)
                continue
            fields = line.split("\t")
            fields.extend([""] * len(col_names))
            new_lines.append("\t".join(fields))

        return new_lines

    # Skip add_column and final_filter if chunked processing handled it
    if not chunked_processing_used:
        if args.add_column:
            buffer = add_named_columns(buffer, args.add_column)

        # Apply final filter if provided
        if cfg.get("final_filter") and len(buffer) > 1:
            logger.info("Applying final filter to results")

            # Convert buffer to DataFrame
            import io

            buffer_str = "\n".join(buffer)
            df = pd.read_csv(io.StringIO(buffer_str), sep="\t", dtype=str, keep_default_na=False)

            # Apply the filter
            from .filters import filter_dataframe_with_query

            df = filter_dataframe_with_query(df, cfg["final_filter"])

            # Convert back to buffer
            output_str = io.StringIO()
            df.to_csv(output_str, sep="\t", index=False, na_rep="")
            buffer = output_str.getvalue().rstrip("\n").split("\n")

            if len(buffer) <= 1:
                logger.warning("No variants passed the final filter")

    if final_to_stdout:
        if buffer:
            sys.stdout.write("\n".join(buffer) + "\n")
        final_out_path = None
    else:
        final_file = final_output
        with open(final_file, "w", encoding="utf-8") as out:
            for line in buffer:
                out.write(line + "\n")
        final_out_path = final_file

    # Gene burden if requested
    if cfg.get("perform_gene_burden", False):
        gene_burden_cfg = cfg.copy()
        gene_burden_cfg["perform_gene_burden"] = True
        line_count = 0
        with open(final_tsv, "r", encoding="utf-8") as inp, open(
            gene_burden_tsv, "w", encoding="utf-8"
        ) as out:
            for line in analyze_variants(inp, gene_burden_cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count <= 1:  # Only header or nothing
            logger.warning("Gene burden requested but no variant data produced.")
        else:
            logger.debug(f"Gene burden analysis complete: {gene_burden_tsv}")

    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Run ended at {end_time.isoformat()}, duration: {duration} seconds")

    # Metadata
    metadata_file = os.path.join(args.output_dir, f"{base_name}.metadata.tsv")
    with open(metadata_file, "w", encoding="utf-8") as mf:
        mf.write("Parameter\tValue\n")

        def meta_write(param, val):
            p = sanitize_metadata_field(param)
            v = sanitize_metadata_field(val)
            mf.write(f"{p}\t{v}\n")

        meta_write("Tool", "variantcentrifuge")
        meta_write("Version", cfg.get("version", "N/A"))
        meta_write("Run_start_time", start_time.isoformat())
        meta_write("Run_end_time", end_time.isoformat())
        meta_write("Run_duration_seconds", str(duration))
        meta_write("Date", datetime.datetime.now().isoformat())
        meta_write("Command_line", " ".join([sanitize_metadata_field(x) for x in sys.argv]))

        for k, v in cfg.items():
            meta_write(f"config.{k}", str(v))

        # Versions
        snpeff_ver = get_tool_version("snpEff")
        bcftools_ver = get_tool_version("bcftools")
        snpsift_ver = get_tool_version("SnpSift")
        bedtools_ver = get_tool_version("bedtools")

        meta_write("tool.snpeff_version", snpeff_ver)
        meta_write("tool.bcftools_version", bcftools_ver)
        meta_write("tool.snpsift_version", snpsift_ver)
        meta_write("tool.bedtools_version", bedtools_ver)

    # Define report directory consistently
    report_dir = os.path.join(args.output_dir, "report")
    os.makedirs(report_dir, exist_ok=True)  # Ensure report_dir exists

    # Phase 1: Excel conversion (without finalization) if requested
    xlsx_file = None
    if args.xlsx and final_out_path and not final_to_stdout:
        if not os.path.exists(final_out_path) or os.path.getsize(final_out_path) == 0:
            logger.warning("Final output file is empty. Cannot convert to Excel.")
        else:
            xlsx_file = convert_to_excel(final_out_path, cfg)
            append_tsv_as_sheet(xlsx_file, metadata_file, sheet_name="Metadata")
            if (
                not cfg["no_stats"]
                and cfg.get("stats_output_file")
                and os.path.exists(cfg["stats_output_file"])
            ):
                if os.path.getsize(cfg["stats_output_file"]) > 0:
                    append_tsv_as_sheet(
                        xlsx_file, cfg["stats_output_file"], sheet_name="Statistics"
                    )
                else:
                    logger.warning("Stats file is empty, skipping Statistics sheet.")
            if (
                cfg.get("perform_gene_burden", False)
                and os.path.exists(gene_burden_tsv)
                and os.path.getsize(gene_burden_tsv) > 0
            ):
                append_tsv_as_sheet(xlsx_file, gene_burden_tsv, sheet_name="Gene Burden")

    # Phase 2: Generate IGV reports and map if IGV is enabled
    # This must happen before finalize_excel_file and produce_report_json
    igv_enabled = cfg.get("igv_enabled", False)
    # Check if the file exists and has content (at least header)
    if igv_enabled and final_out_path and os.path.exists(final_out_path):
        # MODIFIED: Start of local IGV FASTA feature
        bam_map_file = cfg.get("bam_mapping_file")
        igv_reference_genome = cfg.get("igv_reference")
        igv_fasta_file = cfg.get("igv_fasta")
        igv_ideogram_file = cfg.get("igv_ideogram")

        # Validate the required parameters and files
        if not bam_map_file:
            logger.error("For IGV integration, --bam-mapping-file must be provided.")
            sys.exit(1)

        if not igv_fasta_file and not igv_reference_genome:
            logger.error(
                "For IGV integration, either --igv-reference or --igv-fasta must be provided."
            )
            sys.exit(1)

        # Validate local FASTA files if provided
        if igv_fasta_file or igv_ideogram_file:
            from .validators import validate_igv_files

            # Get the FASTA index file if provided
            # FASTA index now follows standard naming convention (FASTA file + .fai)
            validate_igv_files(igv_fasta_file, igv_ideogram_file)

        if igv_fasta_file and igv_reference_genome:
            logger.warning(
                "Both local FASTA file and genome reference ID provided. Local FASTA file will take precedence."
            )
        # MODIFIED: End of local IGV FASTA feature

        from .generate_igv_report import generate_igv_report  # Ensure import is present

        # The output_dir for generate_igv_report should be where igv_reports_map.json is expected
        # which is output_dir/report, and report_path in map are relative to output_dir/report
        # MODIFIED: Start of local IGV FASTA feature
        generate_igv_report(
            variants_tsv=final_out_path,  # final_out_path is the path to the main results TSV
            output_dir=report_dir,  # This ensures map is in report/igv/
            bam_mapping_file=bam_map_file,
            igv_reference=igv_reference_genome,
            integrate_into_main=True,  # Always integrate if IGV is enabled
            igv_fasta=igv_fasta_file,
            igv_ideogram=igv_ideogram_file,
            # MODIFIED: Pass filename shortening parameters
            igv_max_allele_len_filename=cfg.get("igv_max_allele_len_filename", 10),
            igv_hash_len_filename=cfg.get("igv_hash_len_filename", 6),
            igv_max_variant_part_filename=cfg.get("igv_max_variant_part_filename", 50),
            # MODIFIED: Start of IGV flanking feature - pass configurable flanking value
            igv_flanking=cfg.get("igv_flanking", 50),
            # MODIFIED: End of IGV flanking feature
        )
        # MODIFIED: End of local IGV FASTA feature
        logger.info("IGV reports and mapping file generated.")

    # Phase 3: Produce HTML report if requested - after IGV reports are generated
    if args.html_report and final_out_path and os.path.exists(final_out_path):
        # produce_report_json will look for the IGV map created in Phase 2
        produce_report_json(final_out_path, args.output_dir)

        from .generate_html_report import generate_html_report  # Ensure import

        variants_json = os.path.join(report_dir, "variants.json")
        summary_json = os.path.join(report_dir, "summary.json")
        generate_html_report(
            variants_json=variants_json,
            summary_json=summary_json,
            output_dir=report_dir,  # HTML report itself is in report_dir
            cfg=cfg,  # Pass the main configuration dictionary
        )
        logger.info("HTML report generated successfully.")

    # Phase 4: Finalize Excel file after IGV reports are generated
    if xlsx_file:
        # Now finalize the Excel file with awareness of IGV integration
        finalize_excel_file(xlsx_file, cfg)  # cfg contains igv_enabled
        logger.info("Excel file finalized with IGV links (if enabled).")

    # MODIFIED: Start of intermediate cleanup feature
    # Check if intermediate files should be kept (from CLI args or config)
    keep_intermediates = args.keep_intermediates or cfg.get("keep_intermediates", False)

    # Clean up intermediate files if not keeping them
    if not keep_intermediates:
        logger.info("Cleaning up intermediate files...")
        intermediates = [
            variants_file,
            splitted_before_file,
            splitted_after_file,
            filtered_file,
            extracted_tsv,
        ]
        if not args.no_replacement and gt_present:
            intermediates.append(genotype_replaced_tsv)
        if use_phenotypes:
            intermediates.append(phenotype_added_tsv)

        # Safety check: never delete the final output
        if final_out_path and final_out_path in intermediates:
            intermediates.remove(final_out_path)

        # Clean up intermediate directory safely
        intermediate_dir = os.path.join(args.output_dir, "intermediate")
        files_deleted = 0

        for f in intermediates:
            if f and os.path.exists(f):
                try:
                    os.remove(f)
                    files_deleted += 1
                except (IOError, OSError) as e:
                    logger.warning(f"Could not remove intermediate file {f}: {str(e)}")

        logger.info(f"Deleted {files_deleted} intermediate files.")
    else:
        logger.info("Keeping intermediate files (use --keep-intermediates=False to delete them)")
    # MODIFIED: End of intermediate cleanup feature

    logger.info(
        f"Processing completed. Output saved to " f"{'stdout' if final_to_stdout else final_output}"
    )
    
    # Archive results if requested
    if args.archive_results:
        logger.info("Creating compressed archive of results directory...")
        archive_path = archive_results_directory(args.output_dir, base_name)
        if archive_path:
            logger.info(f"Results archived to: {archive_path}")
        else:
            logger.warning("Failed to create results archive, but pipeline completed successfully")

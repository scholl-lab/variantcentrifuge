# File: variantcentrifuge/helpers.py
# Location: variantcentrifuge/variantcentrifuge/helpers.py

"""
Helper functions for analyze_variants and related processes.

Provides:
- Case/control assignment
- Phenotype map building
- Genotype parsing and allele count conversion
- Sample and phenotype classification logic
"""

import logging
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, cast

import pandas as pd

from .phenotype_helpers import (
    assign_case_control_counts,
    build_sample_phenotype_map,
    compute_phenotype_based_case_control_assignment,
    determine_case_control_sets,
    extract_sample_and_genotype,
    genotype_to_allele_count,
)

logger = logging.getLogger("variantcentrifuge")

__all__ = [
    "annotate_variants_with_gene_lists",
    "assign_case_control_counts",
    "build_sample_phenotype_map",
    "check_file",
    "compute_phenotype_based_case_control_assignment",
    "determine_case_control_sets",
    "dump_df_to_xlsx",
    "extract_gencode_id",
    "extract_sample_and_genotype",
    "genotype_to_allele_count",
    "get_vcf_names",
    "get_vcf_regions",
    "get_vcf_samples",
    "get_vcf_size",
    "load_gene_list",
    "match_igv_link_columns",
    "read_sequencing_manifest",
]


def check_file(file_path: str, exit_on_error: bool = True) -> bool:
    """
    Check if a file exists and is readable.

    Parameters
    ----------
    file_path : str
        Path to the file to check
    exit_on_error : bool, optional
        Whether to exit the program if the file does not exist, by default True

    Returns
    -------
    bool
        True if the file exists and is readable, False otherwise

    Raises
    ------
    SystemExit
        If exit_on_error is True and the file does not exist
    """
    if not os.path.exists(file_path):
        if exit_on_error:
            logger.error(f"File not found: {file_path}")
            sys.exit(1)
        return False

    if not os.path.isfile(file_path):
        if exit_on_error:
            logger.error(f"Not a file: {file_path}")
            sys.exit(1)
        return False

    return True


def load_gene_list(file_path: str) -> set[str]:
    """
    Load genes from a file (one gene per line) into a set of uppercase gene names.

    Parameters
    ----------
    file_path : str
        Path to a file containing gene names, one per line

    Returns
    -------
    Set[str]
        A set of gene names in uppercase for case-insensitive matching
    """
    genes: set[str] = set()
    if not os.path.exists(file_path):
        logger.warning(f"Gene list file not found: {file_path}. Skipping this list.")
        return genes

    try:
        with open(file_path, encoding="utf-8") as f:
            for line_num, line_content in enumerate(f, 1):
                gene = line_content.strip()
                if gene:  # Avoid empty lines
                    if (
                        not gene.isalnum() and "_" not in gene and "-" not in gene
                    ):  # Basic check for valid gene name characters
                        logger.debug(
                            f"Potentially non-standard character in gene '{gene}' from file {file_path} at line {line_num}. Including anyway."
                        )
                    genes.add(gene.upper())  # Store as uppercase for case-insensitive matching
    except Exception as e:
        logger.error(f"Error reading gene list file {file_path}: {e}")

    if not genes:
        logger.warning(f"Gene list file {file_path} was empty or no valid genes found.")
    else:
        logger.debug(f"Loaded {len(genes)} genes from {file_path}")

    return genes


def _sanitize_column_name(file_path: str) -> str:
    """
    Convert a file path to a valid TSV column name.

    Extract the filename and remove invalid characters.

    Parameters
    ----------
    file_path : str
        Path to a file

    Returns
    -------
    str
        A sanitized column name derived from the file's basename
    """
    base_name = os.path.basename(file_path)
    name_part, _ext = os.path.splitext(base_name)

    # Iteratively remove common extensions
    common_extensions = [".txt", ".list", ".tsv", ".csv", ".genes", ".gene"]
    original_name_part = name_part
    while True:
        found_ext = False
        for common_ext in common_extensions:
            if name_part.lower().endswith(common_ext):
                name_part = name_part[: -len(common_ext)]
                found_ext = True
                break
        if not found_ext:
            break

    if not name_part:  # If all was extension
        name_part = original_name_part  # Fallback to name before extension stripping

    # Replace non-alphanumeric characters (except underscore) with underscore
    sanitized_name = re.sub(r"[^\w_]", "_", name_part)
    # Remove leading/trailing underscores and collapse multiple underscores
    sanitized_name = re.sub(r"_+", "_", sanitized_name).strip("_")

    if not sanitized_name:  # If name becomes empty after sanitization (e.g. "---.txt")
        sanitized_name = "unnamed_list"

    # Ensure it doesn't start with a number or is purely numeric
    if not sanitized_name or sanitized_name[0].isdigit() or sanitized_name.isdigit():
        sanitized_name = "list_" + sanitized_name

    return sanitized_name


def annotate_variants_with_gene_lists(lines: list[str], gene_list_files: list[str]) -> list[str]:
    """
    Add new columns to variant data lines indicating membership in provided gene lists.

    For each gene list file, a new column is added to the TSV file. The column will contain 'yes' if
    any of the genes in the GENE column (which may be comma-separated) is present in the gene list,
    and 'no' otherwise.

    Parameters
    ----------
    lines : List[str]
        Lines of the TSV file to annotate
    gene_list_files : List[str]
        List of file paths containing gene names (one per line)

    Returns
    -------
    List[str]
        The input lines with additional columns for each gene list

    Notes
    -----
    - If the GENE column is missing, an error is logged and the input lines are returned unchanged
    - Multiple genes in the GENE column can be separated by commas, semicolons, or spaces
    - If a gene list file fails to load, it is skipped with an error message
    - Column names are sanitized from the file basename, duplicates are suffixed with numbers
    """
    if not gene_list_files or not lines or len(lines) < 2:
        return lines

    loaded_gene_sets: dict[str, set[str]] = {}
    ordered_new_column_names: list[str] = []

    for file_path in gene_list_files:
        genes = load_gene_list(file_path)
        if genes:  # Only process if genes were successfully loaded
            base_col_name = _sanitize_column_name(file_path)
            final_col_name = base_col_name
            suffix = 1
            # Ensure unique column names if multiple files sanitize to the same name
            while final_col_name in loaded_gene_sets:
                final_col_name = f"{base_col_name}_{suffix}"
                suffix += 1

            loaded_gene_sets[final_col_name] = genes
            ordered_new_column_names.append(final_col_name)
            logger.info(
                f"Loaded {len(genes)} genes for annotation column '{final_col_name}' from {file_path}"
            )

    if not loaded_gene_sets:
        logger.warning("No valid gene lists were loaded; no annotation columns will be added.")
        return lines

    # Process header
    header_line = lines[0].rstrip("\n")
    header_parts = header_line.split("\t")

    try:
        gene_col_idx = header_parts.index("GENE")
    except ValueError:
        logger.error(
            "CRITICAL: 'GENE' column not found in TSV header. Cannot perform gene list annotation."
        )
        # Return lines unmodified to prevent pipeline breakage if GENE column is unexpectedly missing
        return lines

    new_header_string = "\t".join(header_parts + ordered_new_column_names)
    modified_lines = [new_header_string]

    # Process data lines
    for data_line_str in lines[1:]:
        data_line_str = data_line_str.rstrip("\n")
        if not data_line_str.strip():  # Preserve empty lines if any
            modified_lines.append(data_line_str)
            continue

        fields = data_line_str.split("\t")

        # Ensure fields list is long enough for gene_col_idx
        if gene_col_idx >= len(fields):
            logger.warning(
                f"Line has too few fields to access GENE column (index {gene_col_idx}): {data_line_str}"
            )
            # Append "no" for all new columns for this malformed line
            fields.extend(["" for _ in range(gene_col_idx - len(fields) + 1)])
            fields.extend(["no" for _ in ordered_new_column_names])
            modified_lines.append("\t".join(fields))
            continue

        variant_gene_field_value = fields[gene_col_idx].strip().upper()  # Match in uppercase

        # Handle if GENE column contains multiple genes (e.g., comma-separated or other delimiters)
        if variant_gene_field_value:
            variant_genes_in_record = {
                g.strip() for g in re.split(r"[,; ]+", variant_gene_field_value) if g.strip()
            }
        else:
            variant_genes_in_record = set()

        for col_name in ordered_new_column_names:
            gene_set_for_current_list = loaded_gene_sets[col_name]
            is_member = "no"
            # Check for intersection between variant genes and gene list
            if variant_genes_in_record and not variant_genes_in_record.isdisjoint(
                gene_set_for_current_list
            ):
                is_member = "yes"
            fields.append(is_member)

        modified_lines.append("\t".join(fields))

    logger.info(
        f"Annotation with gene lists complete. Added {len(ordered_new_column_names)} columns."
    )
    return modified_lines


def dump_df_to_xlsx(
    df: pd.DataFrame, output_path: str | Path, sheet_name: str = "Sheet1", index: bool = False
) -> None:
    """
    Export a pandas DataFrame to an Excel (xlsx) file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to export
    output_path : str or Path
        Path where the Excel file will be saved
    sheet_name : str, optional
        Name of the worksheet in the Excel file, by default "Sheet1"
    index : bool, optional
        Whether to include the DataFrame's index in the Excel file, by default False

    Returns
    -------
    None
        The function writes the DataFrame to an Excel file and logs the action
    """
    try:
        output_path = str(output_path)  # Convert Path to string if needed
        df.to_excel(output_path, sheet_name=sheet_name, index=index)
        logger.info(f"DataFrame exported to Excel file: {output_path}")
    except Exception as e:
        logger.error(f"Error exporting DataFrame to Excel: {e}")
        raise


def extract_gencode_id(gene_name: str) -> str:
    """
    Extract the GENCODE ID (ENSG) from a gene name string if present.

    Parameters
    ----------
    gene_name : str
        Gene name, possibly containing an ENSG ID

    Returns
    -------
    str
        The ENSG ID if found, otherwise an empty string
    """
    match = re.search(r"(ENSG\d+)", gene_name)
    return match.group(1) if match else ""


def get_vcf_names(vcf_file: str) -> list[str]:
    """
    Get sample names from a VCF file header.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        List of sample names from the VCF file
    """
    try:
        from subprocess import PIPE, Popen

        cmd = ["bcftools", "query", "-l", vcf_file]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error getting sample names from VCF: {stderr}")
            raise RuntimeError(f"bcftools query -l failed for {vcf_file}: {stderr.strip()}")

        return [line.strip() for line in stdout.split("\n") if line.strip()]
    except RuntimeError:
        raise
    except Exception as e:
        logger.error(f"Error running bcftools query: {e}")
        raise RuntimeError(f"Error running bcftools query on {vcf_file}: {e}") from e


def get_vcf_regions(vcf_file: str) -> list[str]:
    """
    Get the chromosomal regions covered in a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        List of regions in format "chr:start-end"
    """
    try:
        from subprocess import PIPE, Popen

        cmd = ["bcftools", "query", "-f", "%CHROM:%POS\n", vcf_file]
        process = Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error getting regions from VCF: {stderr}")
            return []

        positions = [line.strip() for line in stdout.split("\n") if line.strip()]
        if not positions:
            return []

        # Group by chromosome and find min/max positions
        by_chrom = defaultdict(list)
        for pos in positions:
            if ":" in pos:
                chrom, position = pos.split(":", 1)
                by_chrom[chrom].append(int(position))

        return [
            f"{chrom}:{min(positions)}-{max(positions)}" for chrom, positions in by_chrom.items()
        ]
    except Exception as e:
        logger.error(f"Error analyzing VCF regions: {e}")
        return []


def get_vcf_samples(vcf_file: str) -> list[str]:
    """
    Get an ordered list of sample names from a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    List[str]
        Ordered list of sample names (preserves VCF header order)
    """
    return get_vcf_names(vcf_file)


def get_vcf_size(vcf_file: str) -> int:
    """
    Get the number of variants in a VCF file.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file

    Returns
    -------
    int
        Number of variants in the VCF file
    """
    try:
        from subprocess import PIPE, Popen

        # Use a proper subprocess pipe chain instead of shell=True to
        # prevent command injection via crafted file paths.
        bcf_proc = Popen(
            ["bcftools", "view", "-H", vcf_file],
            stdout=PIPE,
            stderr=PIPE,
        )
        wc_proc = Popen(
            ["wc", "-l"],
            stdin=bcf_proc.stdout,
            stdout=PIPE,
            stderr=PIPE,
            universal_newlines=True,
        )
        # Allow bcf_proc to receive SIGPIPE if wc_proc exits early
        if bcf_proc.stdout:
            bcf_proc.stdout.close()

        stdout, _stderr = wc_proc.communicate()
        bcf_proc.wait()

        if bcf_proc.returncode != 0:
            bcf_stderr = bcf_proc.stderr.read().decode() if bcf_proc.stderr else ""
            logger.error(f"Error getting VCF size: {bcf_stderr}")
            return 0

        try:
            return int(stdout.strip())
        except ValueError:
            logger.error(f"Could not parse VCF size output: {stdout}")
            return 0
    except Exception as e:
        logger.error(f"Error counting variants in VCF: {e}")
        return 0


def match_igv_link_columns(header_fields: list[str]) -> dict[str, int]:
    """
    Find columns in a TSV header that should contain IGV links.

    Parameters
    ----------
    header_fields : List[str]
        List of header field names

    Returns
    -------
    Dict[str, int]
        Mapping of column name to column index for IGV link columns
    """
    # Columns that might contain coordinates for IGV links
    igv_column_patterns = ["CHROM", "POS", "START", "END", "POSITION", "LOCATION"]

    result = {}
    for i, field in enumerate(header_fields):
        if any(pattern in field.upper() for pattern in igv_column_patterns):
            result[field] = i

    return result


def read_sequencing_manifest(file_path: str) -> dict[str, dict[str, str]]:
    """
    Read a sequencing manifest file containing sample metadata.

    Expected format: TSV with header row and sample IDs in first column.

    Parameters
    ----------
    file_path : str
        Path to the manifest file

    Returns
    -------
    Dict[str, Dict[str, str]]
        Mapping of sample IDs to metadata dictionaries
    """
    result: dict[str, dict[str, Any]] = {}
    if not os.path.exists(file_path):
        logger.warning(f"Manifest file not found: {file_path}")
        return result

    try:
        df = pd.read_csv(file_path, sep="\t", dtype=str)
        if df.empty:
            logger.warning(f"Empty manifest file: {file_path}")
            return result

        # Assuming first column contains sample IDs
        id_col = df.columns[0]

        # Convert each row to a dictionary
        for row in df.itertuples(index=True):
            sample_id = getattr(row, id_col, None)
            if sample_id and not pd.isna(sample_id):
                # Convert namedtuple to dict for metadata storage
                result[str(sample_id)] = cast(dict[str, Any], df.loc[row.Index].to_dict())

        logger.info(f"Loaded metadata for {len(result)} samples from {file_path}")
        return result
    except Exception as e:
        logger.error(f"Error reading sequencing manifest: {e}")
        return result

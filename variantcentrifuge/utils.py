# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

Provides helper functions for logging, running commands, checking tool availability,
and retrieving tool versions.
"""

import gzip
import hashlib
import logging
import math
import os
import re
import shutil
import subprocess
from typing import List, Optional, Union

logger = logging.getLogger("variantcentrifuge")


def check_external_tools(tools: List[str]) -> bool:
    """
    Check if external tools are available in PATH.

    Parameters
    ----------
    tools : List[str]
        List of tool names to check for availability

    Returns
    -------
    bool
        True if all tools are available, False otherwise
    """
    for tool in tools:
        if not shutil.which(tool):
            logger.error(f"Required tool not found in PATH: {tool}")
            return False
        logger.debug(f"Found tool in PATH: {tool}")
    return True


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


def run_command(cmd: list, output_file: Optional[str] = None) -> str:
    """
    Run a shell command and write stdout to output_file if provided, else return stdout.

    Parameters
    ----------
    cmd : list of str
        Command and its arguments.
    output_file : str, optional
        Path to a file where stdout should be written. If None,
        returns stdout as a string.

    Returns
    -------
    str
        If output_file is None, returns the command stdout as a string.
        If output_file is provided, returns output_file after completion.

    Raises
    ------
    subprocess.CalledProcessError
        If the command returns a non-zero exit code.
    """
    logger.debug("Running command: %s", " ".join(cmd))
    if output_file:
        with open(output_file, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error("Command failed: %s\nError: %s", " ".join(cmd), result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
    else:
        logger.debug("Command completed successfully.")
        if output_file:
            return output_file
        else:
            return result.stdout


def normalize_vcf_headers(lines: List[str]) -> List[str]:
    """
    Normalize header lines from tools like SnpEff and SnpSift.

    By:
    1. Removing known prefixes (e.g., "ANN[*].", "ANN[0].")
    2. Converting indexed genotype fields from format GEN[index].FIELD to FIELD_index
       (e.g., "GEN[0].AF" -> "AF_0", "GEN[1].DP" -> "DP_1")

    Parameters
    ----------
    lines : List[str]
        A list of lines (e.g., lines from a file) whose first line may contain
        SnpEff/SnpSift-generated prefixes in column headers.

    Returns
    -------
    List[str]
        The updated list of lines where the first line has had matching prefixes
        removed or replaced and indexed fields normalized.
    """
    if not lines:
        return lines

    # Only modify the first line to match the original behavior
    header = lines[0]

    # First apply regex to handle indexed genotype fields (GEN[<index>].<FIELD> -> <FIELD>_<index>)
    # This must be done before other replacements to avoid conflicts
    pattern = r"GEN\[(\d+)\]\.([A-Za-z0-9_]+)"
    header = re.sub(pattern, r"\2_\1", header)
    # Now apply the standard transformations for non-indexed fields
    header = (
        header.replace("ANN[*].", "")
        .replace("ANN[0].", "")
        .replace("GEN[*].", "")
        .replace("GEN[0].", "")
        .replace("NMD[*].", "NMD_")
        .replace("NMD[0].", "NMD_")
        .replace("LOF[*].", "LOF_")
        .replace("LOF[0].", "LOF_")
    )

    lines[0] = header
    return lines


# Keep the original function name as an alias for backward compatibility


def get_tool_version(tool_name: str) -> str:
    """
    Retrieve the version of a given tool.

    Supported tools:

    - snpEff
    - bcftools
    - SnpSift
    - bedtools

    Parameters
    ----------
    tool_name : str
        Name of the tool to retrieve version for.

    Returns
    -------
    str
        Version string or 'N/A' if not found or cannot be retrieved.
    """

    def parse_snpeff(stdout, stderr):
        for line in stdout.splitlines():
            if line.strip():
                return line.strip()
        return "N/A"

    def parse_bcftools(stdout, stderr):
        for line in stdout.splitlines():
            if "bcftools" in line.lower():
                return line.strip()
        return "N/A"

    def parse_snpsift(stdout, stderr):
        for line in stdout.splitlines():
            if line.startswith("SnpSift version"):
                return line.strip()
        for line in stderr.splitlines():
            if line.startswith("SnpSift version"):
                return line.strip()
        return "N/A"

    def parse_bedtools(stdout, stderr):
        for line in stdout.splitlines():
            if "bedtools" in line.lower():
                return line.strip()
        return "N/A"

    tool_map = {
        "snpEff": {"command": ["snpEff", "-version"], "parse_func": parse_snpeff},
        "bcftools": {
            "command": ["bcftools", "--version"],
            "parse_func": parse_bcftools,
        },
        "SnpSift": {"command": ["SnpSift", "annotate"], "parse_func": parse_snpsift},
        "bedtools": {
            "command": ["bedtools", "--version"],
            "parse_func": parse_bedtools,
        },
    }

    if tool_name not in tool_map:
        logger.warning("No version retrieval logic for %s. Returning 'N/A'.", tool_name)
        return "N/A"

    cmd = tool_map[tool_name]["command"]
    parse_func = tool_map[tool_name]["parse_func"]

    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False
        )
        version = parse_func(result.stdout, result.stderr)
        if version == "N/A":
            logger.warning("Could not parse version for %s. Returning 'N/A'.", tool_name)
        return version
    except Exception as e:
        logger.warning("Failed to retrieve version for %s: %s", tool_name, e)
        return "N/A"


def sanitize_metadata_field(value: str) -> str:
    """
    Sanitize a metadata field by removing tabs and newlines, replacing with spaces for TSV.

    Parameters
    ----------
    value : str
        The value to sanitize.

    Returns
    -------
    str
        Sanitized string value with no tabs or newlines.
    """
    # Handle non-string values by converting to string first
    if not isinstance(value, str):
        value = str(value)
    return value.replace("\t", " ").replace("\n", " ").strip()


def ensure_fields_in_extract(base_fields_str: str, extra_fields: List[str]) -> str:
    """
    Ensure each item in extra_fields is present in the space-delimited base_fields_str.

    Notes
    -----
    We no longer normalize extra_fields here, so that raw columns like "GEN[*].DP"
    remain unmodified.
    """
    if not base_fields_str:
        base_list = []
    else:
        base_list = base_fields_str.split()

    # Just deduplicate:
    for raw_field in extra_fields:
        if raw_field not in base_list:
            base_list.append(raw_field)

    return " ".join(base_list)


# MODIFIED: Start of IGV filename shortening feature
def generate_igv_safe_filename_base(
    sample_id: str,
    chrom: str,
    pos: Union[str, int],
    ref: str,
    alt: str,
    max_allele_len: int = 10,
    hash_len: int = 6,
    max_variant_part_len: int = 50,
) -> str:
    """
    Generate a safe, shortened filename base for IGV reports to prevent "File name too long" errors.

    This function handles long REF/ALT alleles by truncating and appending a hash of the original
    allele to maintain uniqueness. The function returns a base filename (without extension) that is
    filesystem-safe and should avoid "File name too long" errors.

    Parameters
    ----------
    sample_id : str
        Sample identifier
    chrom : str
        Chromosome name/identifier
    pos : str or int
        Genomic position
    ref : str
        Reference allele
    alt : str
        Alternate allele
    max_allele_len : int, default=10
        Maximum length for each allele in the filename
    hash_len : int, default=6
        Length of hash to append when truncating an allele
    max_variant_part_len : int, default=50
        Maximum length for the variant part of the filename (chr_pos_ref_alt)

    Returns
    -------
    str
        A safe, shortened filename base
    """
    # Convert position to string if it's an integer
    pos = str(pos)

    # Sanitize inputs (remove characters that might cause issues in filenames)
    sample_id = re.sub(r"[^\w.-]", "_", sample_id)
    chrom = re.sub(r"[^\w.-]", "_", chrom)
    pos = re.sub(r"[^\w.-]", "_", pos)

    # Special case for ref with characters like '[C/G]'
    if any(c in ref for c in ["[", "]", "/", ":", ";"]):
        sanitized_ref = re.sub(r"[^\w.-]", "_", ref)
    else:
        sanitized_ref = ref

    # Special case for alt with special characters
    if any(c in alt for c in ["[", "]", "/", ":", ";"]):
        sanitized_alt = re.sub(r"[^\w.-]", "_", alt)
    else:
        sanitized_alt = alt

    # Handle long ref allele
    if len(ref) > max_allele_len:
        ref_hash = hashlib.md5(ref.encode()).hexdigest()[:hash_len]
        # CRITICAL: Tests expect exactly 3 chars prefix + underscore + 6 char hash
        short_ref = f"{sanitized_ref[:3]}_{ref_hash}"
    else:
        short_ref = sanitized_ref

    # Handle long alt allele
    if len(alt) > max_allele_len:
        alt_hash = hashlib.md5(alt.encode()).hexdigest()[:hash_len]
        # CRITICAL: Tests expect exactly 3 chars prefix + underscore + 6 char hash
        short_alt = f"{sanitized_alt[:3]}_{alt_hash}"
    else:
        short_alt = sanitized_alt

    # Create the variant part of the filename
    variant_part = f"{chrom}_{pos}_{short_ref}_{short_alt}"

    # Check if the variant part is still too long and further shorten if necessary
    if len(variant_part) > max_variant_part_len:
        # Create a hash of the full variant part
        variant_str = f"{chrom}_{pos}_{ref}_{alt}"
        variant_hash = hashlib.md5(variant_str.encode()).hexdigest()[:hash_len]
        # Determine how many chars to keep from each component to fit within max_variant_part_len
        # Reserve space for separators and hash
        # 3 underscores + 1 for separator before hash
        available_chars = max_variant_part_len - len(variant_hash) - 4
        # 4 components: chrom, pos, short_ref, short_alt
        chars_per_component = available_chars // 4

        # Ensure we have at least 1 character per component
        chars_per_component = max(1, chars_per_component)

        # Shorten each component proportionally
        chrom_short = chrom[:chars_per_component]
        pos_short = pos[:chars_per_component]
        ref_short = short_ref[:chars_per_component]
        alt_short = short_alt[:chars_per_component]

        variant_part = f"{chrom_short}_{pos_short}_{ref_short}_{alt_short}_{variant_hash}"

    # Combine sample id and variant part to create the final filename base
    return f"{sample_id}_{variant_part}"


# MODIFIED: End of IGV filename shortening feature


def split_bed_file(input_bed: str, num_chunks: int, output_dir: str) -> List[str]:
    """
    Split a BED file into a specified number of chunks with roughly equal total base pairs.

    Parameters
    ----------
    input_bed : str
        Path to the sorted input BED file.
    num_chunks : int
        The number of smaller BED files to create.
    output_dir : str
        The directory where the chunked BED files will be saved.

    Returns
    -------
    List[str]
        A list of file paths to the created chunked BED files.
    """
    regions = []
    total_bases = 0
    with open(input_bed, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue

            try:
                start, end = int(parts[1]), int(parts[2])
                size = end - start
                regions.append({"line": line.strip(), "size": size})
                total_bases += size
            except (ValueError, IndexError):
                logger.warning(f"Could not parse line in BED file: {line.strip()}")
                continue

    if not regions:
        logger.warning("Input BED file for splitting is empty or invalid.")
        return []

    if num_chunks > len(regions):
        logger.warning(
            f"Number of chunks ({num_chunks}) is greater than number of regions "
            f"({len(regions)}). Reducing chunks to {len(regions)}."
        )
        num_chunks = len(regions)

    target_bases_per_chunk = math.ceil(total_bases / num_chunks)
    logger.info(
        f"Total bases: {total_bases}, splitting into {num_chunks} chunks of "
        f"~{target_bases_per_chunk} bases each."
    )

    chunk_files = []
    current_chunk_bases = 0
    current_chunk_lines = []
    chunk_index = 0

    for region in regions:
        current_chunk_lines.append(region["line"])
        current_chunk_bases += region["size"]

        # If the current chunk is full and we're not on the last chunk, write it.
        if current_chunk_bases >= target_bases_per_chunk and chunk_index < num_chunks - 1:
            chunk_path = os.path.join(output_dir, f"chunk_{chunk_index}.bed")
            with open(chunk_path, "w", encoding="utf-8") as f_chunk:
                f_chunk.write("\n".join(current_chunk_lines) + "\n")
            chunk_files.append(chunk_path)

            # Reset for the next chunk
            chunk_index += 1
            current_chunk_bases = 0
            current_chunk_lines = []

    # Write the last chunk, which contains the remainder
    if current_chunk_lines:
        chunk_path = os.path.join(output_dir, f"chunk_{chunk_index}.bed")
        with open(chunk_path, "w", encoding="utf-8") as f_chunk:
            f_chunk.write("\n".join(current_chunk_lines) + "\n")
        chunk_files.append(chunk_path)

    logger.info(f"Successfully split BED file into {len(chunk_files)} chunks.")
    return chunk_files


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

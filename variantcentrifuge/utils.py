# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

Provides helper functions for logging, running commands, checking tool availability,
and retrieving tool versions.
"""

import hashlib
import logging
import re
import shutil
import subprocess
import sys
from typing import List, Optional, Union

logger = logging.getLogger("variantcentrifuge")

current_log_level = "INFO"


def set_log_level(level: str) -> None:
    """
    Set the global logging level.

    Parameters
    ----------
    level : str
        Logging level, one of ["DEBUG", "INFO", "WARN", "ERROR"].
    """
    global current_log_level
    current_log_level = level


def run_command(cmd: list, output_file: Optional[str] = None) -> str:
    """
    Run a shell command and write its stdout to output_file if provided,
    else return stdout as a string.

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


def normalize_snpeff_headers(lines: List[str]) -> List[str]:
    """
    Remove known SnpEff prefixes (e.g. "ANN[*].", "GEN[0].", etc.) from the first line
    (typically the header line) of the provided lines.

    Parameters
    ----------
    lines : List[str]
        A list of lines (e.g., lines from a file) whose first line may contain
        SnpEff-generated prefixes in column headers.

    Returns
    -------
    List[str]
        The updated list of lines where the first line has had any matching SnpEff
        prefixes removed or replaced.
    """
    if not lines:
        return lines

    # Only modify the first line to match the original behavior
    header = lines[0]
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


def check_external_tools() -> None:
    """
    Check if required external tools are installed and in the PATH.

    Tools checked:
    - bcftools
    - snpEff
    - SnpSift
    - bedtools

    If any are missing, log an error and exit.

    Raises
    ------
    SystemExit
        If any required tool is missing.
    """
    required_tools = ["bcftools", "snpEff", "SnpSift", "bedtools"]
    missing = [tool for tool in required_tools if shutil.which(tool) is None]

    if missing:
        logger.error(
            "Missing required external tools: %s. Please ensure they are installed and in PATH.",
            ", ".join(missing),
        )
        sys.exit(1)
    else:
        logger.debug("All external tools are available.")


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
    Sanitize a metadata field by removing tabs and newlines,
    replacing them with spaces to ensure TSV compatibility.

    Parameters
    ----------
    value : str
        The value to sanitize.

    Returns
    -------
    str
        Sanitized string value with no tabs or newlines.
    """
    return value.replace("\t", " ").replace("\n", " ").strip()


def ensure_fields_in_extract(base_fields_str: str, extra_fields: List[str]) -> str:
    """
    Ensure each item in extra_fields is present in the space-delimited base_fields_str.

    NOTE: We no longer normalize extra_fields here, so that raw columns like "GEN[*].DP"
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

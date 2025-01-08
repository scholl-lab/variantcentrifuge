# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

Provides helper functions for logging, running commands, checking tool availability,
and retrieving tool versions.
"""

import subprocess
import sys
import shutil
import logging
from typing import Optional, List

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
        header
        .replace("ANN[*].", "")
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
            ", ".join(missing)
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
        "snpEff": {
            "command": ["snpEff", "-version"],
            "parse_func": parse_snpeff
        },
        "bcftools": {
            "command": ["bcftools", "--version"],
            "parse_func": parse_bcftools
        },
        "SnpSift": {
            "command": ["SnpSift", "annotate"],
            "parse_func": parse_snpsift
        },
        "bedtools": {
            "command": ["bedtools", "--version"],
            "parse_func": parse_bedtools
        }
    }

    if tool_name not in tool_map:
        logger.warning("No version retrieval logic for %s. Returning 'N/A'.", tool_name)
        return "N/A"

    cmd = tool_map[tool_name]["command"]
    parse_func = tool_map[tool_name]["parse_func"]

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
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
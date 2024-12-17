# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

Provides helper functions for logging, running commands, and checking tool availability,
and retrieving tool versions.
"""

import subprocess
import sys
import shutil
import logging

logger = logging.getLogger("variantcentrifuge")

current_log_level = "INFO"


def set_log_level(level):
    """
    Set the global logging level.

    Parameters
    ----------
    level : str
        Logging level, one of ["DEBUG", "INFO", "WARN", "ERROR"].
    """
    global current_log_level
    current_log_level = level


def run_command(cmd, output_file=None):
    """
    Run a shell command and write its stdout to output_file if provided,
    else return stdout.

    Parameters
    ----------
    cmd : list
        Command and arguments.
    output_file : str, optional
        Path to a file where stdout should be written. If None,
        return stdout as a string.

    Returns
    -------
    str
        If output_file is None, returns the command stdout as a string.
        If output_file is given, returns output_file after completion.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    if output_file:
        with open(output_file, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f,
                                    stderr=subprocess.PIPE,
                                    text=True)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        logger.error(f"Command failed: {' '.join(cmd)}\nError: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
    else:
        logger.debug("Command completed successfully.")
        if output_file:
            return output_file
        else:
            return result.stdout


def check_external_tools():
    """
    Check if required external tools are installed and in the PATH.

    Tools checked:
    - bcftools
    - snpEff
    - SnpSift
    - bedtools

    If any are missing, an error is logged and the program exits.
    """
    required_tools = ["bcftools", "snpEff", "SnpSift", "bedtools"]
    missing = [tool for tool in required_tools if shutil.which(tool) is None]

    if missing:
        logger.error(f"Missing required external tools: {', '.join(missing)}. "
                     f"Please ensure they are installed and in PATH.")
        sys.exit(1)
    else:
        logger.debug("All external tools are available.")


def get_tool_version(tool_name):
    """
    Retrieve the version of a given tool.

    Parameters
    ----------
    tool_name : str
        Name of the tool to retrieve version for.
        Supported: snpEff, bcftools, SnpSift, bedtools

    Returns
    -------
    str
        Version string or 'N/A' if not found or cannot be retrieved.
    """

    # Define how to retrieve versions for specific tools
    # Each entry: tool_name: { "command": [...], "parse_func": callable }
    # parse_func should accept (stdout, stderr) and return a version string or "N/A"
    def parse_snpeff(stdout, stderr):
        # snpEff -version
        # Usually prints something like: 4.3t (build 2017-11-24)
        # We'll take the first line of stdout if present
        for line in stdout.splitlines():
            if line.strip():
                return line.strip()
        return "N/A"

    def parse_bcftools(stdout, stderr):
        # bcftools --version
        # Usually prints something like: bcftools 1.9-74-g15f4c3c
        for line in stdout.splitlines():
            if "bcftools" in line.lower():
                return line.strip()
        return "N/A"

    def parse_snpsift(stdout, stderr):
        # SnpSift annotate
        # Typically prints help text including "SnpSift version x.x"
        # Check stdout first, then stderr
        for line in stdout.splitlines():
            if line.startswith("SnpSift version"):
                return line.strip()
        for line in stderr.splitlines():
            if line.startswith("SnpSift version"):
                return line.strip()
        return "N/A"

    def parse_bedtools(stdout, stderr):
        # bedtools --version
        # Typically prints something like: bedtools v2.29.2
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
        logger.warning(f"No version retrieval logic for {tool_name}. Returning 'N/A'.")
        return "N/A"

    cmd = tool_map[tool_name]["command"]
    parse_func = tool_map[tool_name]["parse_func"]

    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, check=False)
        version = parse_func(result.stdout, result.stderr)
        if version == "N/A":
            logger.warning(f"Could not parse version for {tool_name}. Returning 'N/A'.")
        return version
    except Exception as e:
        logger.warning(f"Failed to retrieve version for {tool_name}: {e}")
        return "N/A"

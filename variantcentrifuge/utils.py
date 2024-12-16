# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

Provides helper functions for logging, running commands, and checking tool availability.
"""

import subprocess
import sys
import shutil
import tempfile
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
    Run a shell command and write its stdout to output_file if provided, else return stdout.

    Parameters
    ----------
    cmd : list
        Command and arguments.
    output_file : str, optional
        Path to a file where stdout should be written. If None, return stdout as string.

    Returns
    -------
    str
        If output_file is None, returns the command stdout as a string.
        If output_file is given, returns output_file after completion.
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    if output_file:
        with open(output_file, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    else:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
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

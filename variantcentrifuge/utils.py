# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

This module provides helper functions for logging and running external
commands and checking for external tool availability.
"""

import subprocess
import sys
import shutil

def log_message(level, message):
    """
    Log messages to stderr with a given level.

    Parameters
    ----------
    level : str
        Log level (e.g. "INFO", "ERROR").
    message : str
        Message to log.
    """
    print(f"[{level}] {message}", file=sys.stderr)


def run_command_stream(cmd, input_stream=None):
    """
    Run a command and yield its stdout lines.

    Parameters
    ----------
    cmd : list
        Command and arguments as a list.
    input_stream : iterator, optional
        If provided, lines from this iterator will be passed to the
        subprocess stdin.

    Yields
    ------
    str
        Lines from the command's stdout.
    """
    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE if input_stream else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

    if input_stream:
        for line in input_stream:
            proc.stdin.write(line + "\n")
        proc.stdin.close()

    for out_line in proc.stdout:
        yield out_line.rstrip("\n")

    proc.stdout.close()
    proc.wait()

    if proc.returncode != 0:
        err = proc.stderr.read()
        proc.stderr.close()
        raise subprocess.CalledProcessError(
            proc.returncode, cmd, err
        )
    else:
        proc.stderr.close()


def check_external_tools():
    """
    Check if required external tools are installed and in the PATH.

    The tools checked are:
    - bcftools
    - snpEff
    - SnpSift
    - bedtools

    If any are missing, log an error and exit.
    """
    required_tools = ["bcftools", "snpEff", "SnpSift", "bedtools"]
    missing = [tool for tool in required_tools if shutil.which(tool) is None]

    if missing:
        log_message("ERROR", f"Missing required external tools: {', '.join(missing)}. "
                             f"Please ensure they are installed and in PATH.")
        sys.exit(1)

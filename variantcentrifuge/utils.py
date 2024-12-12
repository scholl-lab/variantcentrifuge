# File: variantcentrifuge/utils.py
# Location: variantcentrifuge/variantcentrifuge/utils.py

"""
Utility functions module.

This module provides helper functions for logging and running external
commands.
"""

import subprocess
import sys


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

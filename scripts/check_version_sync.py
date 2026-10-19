#!/usr/bin/env python3
"""Check synchronization of version metadata across repository files."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


def extract_python_version(path: Path) -> str:
    """Extract __version__ from variantcentrifuge/version.py."""
    if not path.exists():
        raise FileNotFoundError(f"Python version file not found: {path}")

    content = path.read_text(encoding="utf-8")
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', content, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find __version__ in {path}")
    return match.group(1).strip()


def extract_cff_version(path: Path) -> str:
    """Extract version from CITATION.cff."""
    if not path.exists():
        raise FileNotFoundError(f"CITATION.cff not found: {path}")

    content = path.read_text(encoding="utf-8")
    match = re.search(r"^version:\s*[\"']?([^\"'\s]+)[\"']?", content, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find version in {path}")
    return match.group(1).strip()


def extract_zenodo_version(path: Path) -> str:
    """Extract version from .zenodo.json."""
    if not path.exists():
        raise FileNotFoundError(f".zenodo.json not found: {path}")

    content = path.read_text(encoding="utf-8")
    data = json.loads(content)
    version = data.get("version")
    if not version:
        raise ValueError(f"Could not find 'version' key in {path}")
    return str(version).strip()


def check_version_sync(repo_root: Path, tag: str | None = None) -> dict[str, str]:
    """
    Verify that version.py, CITATION.cff, and .zenodo.json have matching versions.

    Parameters
    ----------
    repo_root : Path
        Repository root directory.
    tag : str | None
        Optional git tag (e.g. 'v0.17.16' or '0.17.16') to validate against files.

    Returns
    -------
    dict[str, str]
        Mapping from file identifier to extracted version string.

    Raises
    ------
    ValueError
        If versions across files disagree or if tag does not match file versions.
    """
    version_py = repo_root / "variantcentrifuge" / "version.py"
    citation_cff = repo_root / "CITATION.cff"
    zenodo_json = repo_root / ".zenodo.json"

    py_ver = extract_python_version(version_py)
    cff_ver = extract_cff_version(citation_cff)
    zenodo_ver = extract_zenodo_version(zenodo_json)

    versions = {
        "version.py": py_ver,
        "CITATION.cff": cff_ver,
        ".zenodo.json": zenodo_ver,
    }

    unique_versions = set(versions.values())
    if len(unique_versions) > 1:
        details = ", ".join(f"{k}={v}" for k, v in versions.items())
        raise ValueError(f"Version mismatch detected: {details}")

    canonical_version = py_ver

    if tag is not None:
        clean_tag = tag.removeprefix("v")
        if clean_tag != canonical_version:
            raise ValueError(
                f"Tag version '{clean_tag}' does not match canonical version '{canonical_version}'"
            )

    return versions


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Verify version synchronization across version.py, CITATION.cff, and .zenodo.json"
        ),
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Path to repository root (defaults to parent of scripts/)",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default=None,
        help="Optional release tag (e.g. 'v0.17.16') to validate against metadata files",
    )
    return parser.parse_args(args)


def main(args: list[str] | None = None) -> int:
    """Entry point for checking version synchronization."""
    parsed = parse_args(args)
    try:
        versions = check_version_sync(parsed.repo_root, tag=parsed.tag)
        canonical = versions["version.py"]
        print(f"Version synchronization check passed: {canonical}")
        return 0
    except Exception as exc:
        print(f"Error checking version synchronization: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

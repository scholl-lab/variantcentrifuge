#!/usr/bin/env python3
"""Enforce the production Python module line-count budget."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

MAX_LINES = 600


def parse_allowlist(path: Path | None) -> dict[Path, int]:
    """Return resolved source paths and their permitted line-count ceilings."""
    if path is None:
        return {}

    entries: dict[Path, int] = {}
    for number, raw_line in enumerate(path.read_text().splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            raw_path, raw_ceiling = line.rsplit(":", maxsplit=1)
            ceiling = int(raw_ceiling)
        except ValueError as error:
            message = f"Invalid allowlist entry at {path}:{number}: {raw_line!r}"
            raise ValueError(message) from error
        if ceiling <= MAX_LINES:
            message = f"Allowlist ceiling must exceed {MAX_LINES} at {path}:{number}"
            raise ValueError(message)
        entries[Path(raw_path).resolve()] = ceiling
    return entries


def count_lines(path: Path) -> int:
    """Return the physical line count for a UTF-8 Python source file."""
    with path.open(encoding="utf-8") as source_file:
        return sum(1 for _ in source_file)


def is_test_module(path: Path, source_root: Path) -> bool:
    """Return whether a source path is under a directory named ``tests``."""
    return "tests" in path.relative_to(source_root).parts


def find_violations(source_root: Path, allowlist: dict[Path, int]) -> list[str]:
    """Return all files that exceed their applicable line-count ceiling."""
    violations: list[str] = []
    for path in sorted(source_root.rglob("*.py")):
        if is_test_module(path, source_root):
            continue
        lines = count_lines(path)
        ceiling = allowlist.get(path.resolve(), MAX_LINES)
        if lines > ceiling:
            violations.append(f"{path}: {lines} lines exceeds the {ceiling}-line limit")
    return violations


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check production Python modules against the line-count budget."
    )
    parser.add_argument("source_root", type=Path, help="Production source directory")
    parser.add_argument(
        "--allowlist",
        type=Path,
        help="Optional path:maximum_lines baseline allowlist",
    )
    return parser.parse_args()


def main() -> int:
    """Run the line-count check and return a process exit status."""
    args = parse_args()
    source_root = args.source_root.resolve()
    if not source_root.is_dir():
        print(f"Source root is not a directory: {args.source_root}", file=sys.stderr)
        return 2

    try:
        allowlist = parse_allowlist(args.allowlist)
    except (OSError, ValueError) as error:
        print(error, file=sys.stderr)
        return 2

    violations = find_violations(source_root, allowlist)
    if violations:
        print("Python module line-count violations:")
        print("\n".join(violations))
        return 1

    print(f"Line-count check passed for {source_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

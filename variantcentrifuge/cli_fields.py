"""Field file loading utilities for VariantCentrifuge CLI."""

from __future__ import annotations

from pathlib import Path


def load_fields_from_file(fields_file: str | Path) -> str:
    """Read fields to extract from a text file, ignoring blank lines and comments."""
    path = Path(fields_file)
    if not path.is_file():
        raise FileNotFoundError(f"Fields file not found: {fields_file}")
    fields: list[str] = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                fields.extend(p.strip() for p in stripped.replace(",", " ").split() if p.strip())
    return " ".join(fields)

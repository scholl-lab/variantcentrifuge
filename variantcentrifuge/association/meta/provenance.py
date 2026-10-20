# File: variantcentrifuge/association/meta/provenance.py
"""
Annotation provenance for score artifacts (#197).

Recorded at emission from the input VCF's header (``##SnpEffVersion``,
``##SnpEffCmd``, ``##SnpSiftCmd``/other ``*Cmd``/``*Command`` lines) and SHA-256
hashes of the ``--regions-bed`` / ``--gene-file`` contents; compared across strata
by ``meta``. ``None`` means unknown (absent header line or an older artifact) and
only warns; a file hash of ``"none"`` means no such file was used.
"""

from __future__ import annotations

import gzip
import hashlib
import logging
import os
from typing import Any

logger = logging.getLogger("variantcentrifuge")

PROVENANCE_FIELDS = (
    "snpeff_version",
    "snpeff_genome_db",
    "dbnsfp_version",
    "gnomad_source",
    "regions_bed_sha256",
    "gene_file_sha256",
)


def _vcf_meta_lines(path: str) -> list[str]:
    """The ``##`` header lines of a plain or (b)gzipped VCF, or a BCF."""
    with open(path, "rb") as fh:
        is_gz = fh.read(2) == b"\x1f\x8b"
    with gzip.open(path, "rb") if is_gz else open(path, "rb") as fh:
        if fh.read(3) == b"BCF":  # BCF2: magic(5) + uint32 l_text + NUL-terminated header
            fh.read(2)
            size = int.from_bytes(fh.read(4), "little")
            text = fh.read(size).rstrip(b"\0").decode("utf-8", "replace")
            return [ln for ln in text.splitlines() if ln.startswith("##")]
        fh.seek(0)
        lines: list[str] = []
        for raw in fh:
            line = raw.decode("utf-8", "replace").rstrip("\r\n")
            if not line.startswith("##"):
                break
            lines.append(line)
    return lines


def _header_value(lines: list[str], key: str) -> list[str]:
    prefix = f"##{key}="
    return [ln[len(prefix) :].strip().strip('"').strip() for ln in lines if ln.startswith(prefix)]


def _snpeff_genome(cmd: str) -> str | None:
    # SnpEff [cmd] [opts] genome [vcf]: the genome is the last positional before the VCF.
    # ponytail: an option value placed after the genome (rare) would be misread.
    args = [t for t in cmd.split()[1:] if not t.startswith("-")]
    args = [t for t in args if t not in ("eff", "ann")]
    if args and (".vcf" in args[-1] or args[-1] == "-"):
        args = args[:-1]
    return args[-1] if args else None


def _command_files(lines: list[str], needle: str) -> str | None:
    """Basenames of the annotation databases in ``##*Cmd=``/``##*Command=`` lines that
    contain ``needle``: the ``-db`` value, else the first file argument (SnpSift annotate
    and bcftools annotate -a put the database first; the input VCF comes after it).
    """
    names: set[str] = set()
    for ln in lines:
        key, _, value = ln[2:].partition("=")
        if not key.endswith(("Cmd", "Command")):
            continue
        toks = value.strip().strip('"').split()
        if "-db" in toks[:-1]:
            db: str | None = toks[toks.index("-db") + 1]
        else:
            db = next((t for t in toks[1:] if not t.startswith("-") and "." in t), None)
        if db and needle in db.lower():
            names.add(os.path.basename(db))
    return ",".join(sorted(names)) or None


def _file_sha256(path: str | None) -> str:
    if not path:
        return "none"
    with open(path, "rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def collect_provenance(config: dict[str, Any]) -> dict[str, str | None]:
    """Annotation provenance of this run (see module docstring for the conventions)."""
    vcf = config.get("vcf_file")
    lines = _vcf_meta_lines(vcf) if vcf and os.path.isfile(vcf) else []
    versions = _header_value(lines, "SnpEffVersion")
    cmds = _header_value(lines, "SnpEffCmd")
    return {
        "snpeff_version": versions[-1] if versions else None,
        "snpeff_genome_db": _snpeff_genome(cmds[-1]) if cmds else None,
        "dbnsfp_version": _command_files(lines, "dbnsfp"),
        "gnomad_source": _command_files(lines, "gnomad"),
        "regions_bed_sha256": _file_sha256(config.get("regions_bed")),
        "gene_file_sha256": _file_sha256(config.get("gene_file")),
    }


def check_provenance(artifacts: list, allow_mismatch: bool, error: type[Exception]) -> None:
    """Refuse (``error``) or, with ``allow_mismatch``, warn when strata differ; warn on unknown."""
    unknown: list[str] = []
    for name in PROVENANCE_FIELDS:
        values = {a.stratum_id: a.provenance.get(name) for a in artifacts}
        missing = [s for s, v in values.items() if v is None]
        if missing:
            unknown.append(f"{name} ({','.join(missing)})")
        if len({v for v in values.values() if v is not None}) > 1:
            listing = ", ".join(f"{s}={v}" for s, v in values.items() if v is not None)
            msg = f"Strata disagree on annotation provenance '{name}': {listing}."
            if not allow_mismatch:
                raise error(f"{msg} Pass --allow-annotation-mismatch to pool anyway.")
            logger.warning("%s Pooling anyway (--allow-annotation-mismatch).", msg)
    if unknown:
        logger.warning("Annotation provenance unknown, so not checked: %s.", "; ".join(unknown))


__all__ = ["PROVENANCE_FIELDS", "check_provenance", "collect_provenance"]

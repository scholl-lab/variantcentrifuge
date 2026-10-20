"""Fail-fast SnpSift check of ``--filters`` on the VCF header plus one synthetic record (#211)."""

from __future__ import annotations

import logging
import os
import re
import subprocess
import tempfile

logger = logging.getLogger(__name__)


def _synthetic_vcf(header: str) -> str:
    """Header plus one record built only from the header (no real record is read)."""
    lines = header.rstrip("\n").split("\n")
    contig = next(
        (m.group(1) for ln in lines if (m := re.match(r"##contig=<ID=([^,>]+)", ln))), "1"
    )
    samples = lines[-1].split("\t")[9:]
    rec = [contig, "1", ".", "A", "C", ".", ".", "."]
    if samples:
        rec += ["GT"] + ["./."] * len(samples)
    return "\n".join([*lines, "\t".join(rec)]) + "\n"


def validate_snpsift_filter(vcf_file: str, filter_expr: str, snpsift_flags: list[str]) -> None:
    """Raise ValueError with SnpSift's own message if it cannot evaluate ``filter_expr``.

    ponytail: SnpSift resolves fields lazily on the first record, hence the synthetic record;
    a field only reached behind a short-circuited ``&``/``|`` is not caught here.
    """
    header = subprocess.run(
        ["bcftools", "view", "-h", vcf_file], capture_output=True, text=True, check=True
    ).stdout
    fd, path = tempfile.mkstemp(suffix=".vcf", prefix="filter_dry_run_")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(_synthetic_vcf(header))
        res = subprocess.run(
            ["SnpSift", *snpsift_flags, "filter", filter_expr, path],
            capture_output=True,
            text=True,
        )
    finally:
        os.remove(path)
    if res.returncode != 0:
        raise ValueError(
            f"SnpSift rejected the filter expression: {filter_expr}\n"
            f"{res.stderr.strip()[:1500]}\n"
            "Hint: string literals must be quoted, e.g. ANN[0].IMPACT = 'HIGH', and every "
            "field must exist in the VCF header."
        )
    logger.debug("Filter expression accepted by SnpSift dry run: %s", filter_expr)

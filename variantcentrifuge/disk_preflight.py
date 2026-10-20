"""Fail-fast free-space check on the temp and output filesystems (#213)."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

# Measured on tests/fixtures/1000g_cohort: the uncompressed VCF is 15.9x the .vcf.gz and the
# extracted TSV 10.3x; the largest intermediates coexist, so warn below ~26x, rounded up.
WARN_MULTIPLE = 30


def check_free_space(vcf_file: str | os.PathLike, dirs: list[Path]) -> None:
    """Raise OSError if a filesystem has less free space than the input VCF; warn below 30x."""
    if not vcf_file or not os.path.isfile(vcf_file):
        return
    size = os.path.getsize(vcf_file)
    for d in dict.fromkeys(Path(x) for x in dirs):
        free = shutil.disk_usage(d).free
        if free < size:
            raise OSError(
                f"Only {free / 1e9:.2f} GB free in {d} but the input VCF is {size / 1e9:.2f} GB "
                "compressed; point --temp-dir / --output-dir at a larger filesystem."
            )
        if free < WARN_MULTIPLE * size:
            logger.warning(
                "Free space in %s (%.1f GB) is below %dx the compressed input VCF (%.2f GB); "
                "the run may fail with 'No space left on device'. Consider --temp-dir.",
                d,
                free / 1e9,
                WARN_MULTIPLE,
                size / 1e9,
            )

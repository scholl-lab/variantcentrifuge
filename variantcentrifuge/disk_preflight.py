"""Advisory free-space check on the temp and output filesystems (#213)."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

# On tests/fixtures/1000g_cohort the uncompressed VCF is 15.9x the .vcf.gz (the largest single
# intermediate of a whole-VCF run; the extracted TSV adds 10.3x). Gene-panel runs need far less,
# so this is only a warning, and only below the single largest intermediate.
WARN_MULTIPLE = 16


def check_free_space(vcf_file: str | os.PathLike, dirs: list[Path]) -> None:
    """Warn when a filesystem has less than WARN_MULTIPLE x the compressed input VCF free."""
    if not vcf_file or not os.path.isfile(vcf_file):
        return
    size = os.path.getsize(vcf_file)
    for d in dict.fromkeys(Path(x) for x in dirs):
        free = shutil.disk_usage(d).free
        if free < WARN_MULTIPLE * size:
            logger.warning(
                "Free space in %s (%.1f GB) is below %dx the compressed input VCF (%.2f GB); "
                "a whole-VCF run may fail with 'No space left on device'. Consider --temp-dir.",
                d,
                free / 1e9,
                WARN_MULTIPLE,
                size / 1e9,
            )

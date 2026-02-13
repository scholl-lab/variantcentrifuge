"""Pure helper functions for the variantcentrifuge Snakemake workflow.

This module deliberately avoids importing anything from Snakemake so that it
can be unit-tested in a plain ``pytest`` session without Snakemake installed.
"""

from __future__ import annotations

import os
import shlex
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd


def get_samples(samples_df: pd.DataFrame) -> list[str]:
    """Return a sorted list of unique sample names from the sample sheet."""
    return sorted(samples_df["sample"].unique().tolist())


def get_vcf_path(sample: str, vcf_folder: str, samples_df: pd.DataFrame) -> str:
    """Resolve the full VCF path for *sample* using the sample sheet."""
    row = samples_df.loc[sample]
    basename: str = row["vcf_basename"]
    return os.path.join(vcf_folder, f"{basename}.vcf.gz")


def build_vc_flags(vc_config: dict[str, Any], igv_config: dict[str, Any]) -> list[str]:
    """Build variantcentrifuge CLI flags from nested config dicts.

    Returns a list of flag strings that can be joined with ``" ".join()``.
    """
    flags: list[str] = []

    if vc_config.get("add_chr"):
        flags.append("--add-chr")
    if vc_config.get("generate_xlsx"):
        flags.append("--xlsx")
    if vc_config.get("generate_html_report"):
        flags.append("--html-report")

    profile = vc_config.get("field_profile")
    if profile:
        flags.append(f"--field-profile {shlex.quote(str(profile))}")

    for preset in vc_config.get("presets", []):
        flags.append(f"--preset {shlex.quote(str(preset))}")

    gt_fields: list[str] = vc_config.get("append_genotype_fields", [])
    if gt_fields:
        quoted = " ".join(shlex.quote(f) for f in gt_fields)
        flags.append(f"--append-extra-sample-fields {quoted}")

    for gl in vc_config.get("gene_list_files", []):
        flags.append(f"--annotate-gene-list {shlex.quote(gl)}")

    # IGV flags
    if igv_config.get("enabled"):
        flags.append("--igv")
        ref = igv_config.get("reference")
        if ref:
            flags.append(f"--igv-reference {shlex.quote(str(ref))}")
        bam = igv_config.get("bam_mapping_file")
        if bam:
            flags.append(f"--bam-mapping-file {shlex.quote(str(bam))}")
        fasta = igv_config.get("fasta")
        if fasta:
            flags.append(f"--igv-fasta {shlex.quote(str(fasta))}")
        ideogram = igv_config.get("ideogram")
        if ideogram:
            flags.append(f"--igv-ideogram {shlex.quote(str(ideogram))}")
        flags.append(f"--igv-flanking {shlex.quote(str(igv_config.get('flanking', 50)))}")

    return flags

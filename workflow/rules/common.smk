"""Common variables and helper wrappers for the variantcentrifuge workflow."""

import os

from rules.helpers import build_vc_flags, get_samples, get_vcf_path  # noqa: F401

# ── Config shortcuts ─────────────────────────────────────────────────────────
VCF_FOLDER = config["paths"]["vcf_folder"]
OUTPUT_DIR = config["paths"]["output_folder"]
LOG_SUBDIR = config["paths"].get("log_subdir", "logs")
VC_CONFIG = config["variantcentrifuge"]
IGV_CONFIG = config.get("igv", {})
CONTAINER_IMAGE = config.get("container", {}).get("image", "")

SAMPLES = get_samples(samples_df)


def get_final_outputs():
    """Return the list of final output files expected by *rule all*."""
    return expand(
        os.path.join(OUTPUT_DIR, "{sample}", "{sample}.vc_analysis.tsv"),
        sample=SAMPLES,
    )

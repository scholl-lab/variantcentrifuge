"""Helper functions for checkpoint configuration hashing, validation, and reporting."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
from typing import Any

from .config import normalize_annotation_config
from .sample_subset import blocking_reason, sample_set_sha256

NEW_FINGERPRINT_KEYS: tuple[str, ...] = (
    "filter_exception_file",
    "transcript_list",
    "transcript_file",
    "regions_bed",
    "interval_expand",
    "add_chr",
    "tumor_sample_index",
    "normal_sample_index",
)

RELEVANT_CONFIG_KEYS: tuple[str, ...] = (
    "gene_name",
    "gene_file",
    "vcf_file",
    "output_file",
    "preset",
    "filters",
    "fields_to_extract",
    "field_profile",
    "threads",
    "no_replacement",
    "late_filtering",
    "phenotype_file",
    "phenotype_sample_column",
    "phenotype_value_column",
    "case_phenotypes",
    "control_phenotypes",
    "case_samples_file",
    "control_samples_file",
    "case_samples",
    "control_samples",
    "no_sample_subset",
    "reference",
    "sample_sheet",
    "from_tsv",
    "bcftools_prefilter",
    "final_filter",
    "scoring_config_path",
    "ped",
    "inheritance_mode",
    "annotate_bed_files",
    "annotate_gene_lists",
    "annotate_json_genes",
    "filter_exception_file",
    "transcript_list",
    "transcript_file",
    "regions_bed",
    "interval_expand",
    "add_chr",
    "tumor_sample_index",
    "normal_sample_index",
    "gene_set_file",
    "mask_label",
)

CANONICAL_ANNOTATION_KEYS: frozenset[str] = frozenset(
    {
        "annotate_bed_files",
        "annotate_gene_lists",
    }
)

LEGACY_RELEVANT_CONFIG_KEYS: tuple[str, ...] = (
    "gene_name",
    "gene_file",
    "vcf_file",
    "output_file",
    "preset",
    "filters",
    "fields",
    "threads",
    "no_replacement",
    "late_filtering",
    "chunk_size",
    "samples_file",
    "bcftools_prefilter",
    "final_filter",
    "scoring_config_path",
    "ped",
    "inheritance_mode",
    "annotate_bed",
    "annotate_gene_list",
    "annotate_json_genes",
)

ANNOTATION_ALIAS_KEYS: frozenset[str] = frozenset(
    {
        "annotate_bed",
        "annotate_bed_files",
        "annotate_gene_list",
        "annotate_gene_list_files",
        "annotate_gene_lists",
    }
)


# Bump when BED overlap semantics change so old checkpoints with BED annotations are not reused.
BED_OVERLAP_VERSION = 2


def hash_configuration(config: dict[str, Any]) -> str:
    """Create a hash of the configuration for change detection."""
    norm_config = normalize_annotation_config(config.copy())

    relevant_config = {
        k: v
        for k, v in norm_config.items()
        if k in RELEVANT_CONFIG_KEYS
        and v is not None
        and not (k in CANONICAL_ANNOTATION_KEYS and not v)
    }
    relevant_config["__fingerprint_version__"] = 3
    if relevant_config.get("annotate_bed_files"):
        relevant_config["__bed_overlap__"] = BED_OVERLAP_VERSION
    if relevant_config.get("gene_set_file"):  # edited in place must also invalidate resume
        with open(relevant_config["gene_set_file"], "rb") as fh:
            relevant_config["__gene_set_sha256__"] = hashlib.sha256(fh.read()).hexdigest()
    # A blocker flag (--no-stats, --pseudonymize, ...) decides whether samples are subset (#212).
    if config.get("perform_association"):
        relevant_config["__sample_subset_blocked__"] = blocking_reason(config) is not None
    # Edited sample/phenotype files change the extracted sample subset (#212).
    for key in ("case_samples_file", "control_samples_file", "phenotype_file"):
        path = relevant_config.get(key)
        if path and os.path.isfile(path):
            with open(path, "rb") as fh:
                data = fh.read()
            relevant_config[f"__{key}_sha256__"] = (
                sample_set_sha256([ln.strip() for ln in data.decode().splitlines() if ln.strip()])
                if key != "phenotype_file"
                else hashlib.sha256(data).hexdigest()
            )

    config_str = json.dumps(relevant_config, sort_keys=True)
    return hashlib.sha256(config_str.encode()).hexdigest()


def hash_legacy_configuration(config: dict[str, Any]) -> str:
    """Create the pre-normalization checkpoint hash for compatibility."""
    relevant_config = {
        k: v for k, v in config.items() if k in LEGACY_RELEVANT_CONFIG_KEYS and v is not None
    }
    config_str = json.dumps(relevant_config, sort_keys=True)
    return hashlib.sha256(config_str.encode()).hexdigest()


def legacy_annotation_equivalent_configs(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Build legacy-key configs that represent the normalized annotation behavior."""
    normalized = normalize_annotation_config(config.copy())
    base_config = config.copy()
    for key in ANNOTATION_ALIAS_KEYS:
        base_config.pop(key, None)

    bed_files = normalized.get("annotate_bed_files", [])
    gene_lists = normalized.get("annotate_gene_lists", [])

    bed_options = [{"annotate_bed": bed_files}] if bed_files else [{}, {"annotate_bed": []}]
    gene_options = (
        [{"annotate_gene_list": gene_lists}] if gene_lists else [{}, {"annotate_gene_list": []}]
    )

    legacy_configs = []
    for bed_option in bed_options:
        for gene_option in gene_options:
            candidate = base_config.copy()
            candidate.update(bed_option)
            candidate.update(gene_option)
            legacy_configs.append(candidate)
    return legacy_configs


def resume_configuration_hashes(config: dict[str, Any]) -> set[str]:
    """Return current and legacy-annotation-key-equivalent config hashes accepted for resume."""
    hashes = {hash_configuration(config)}
    for candidate in legacy_annotation_equivalent_configs(config):
        hashes.add(hash_configuration(candidate))
    return hashes


def format_state_summary(state: dict[str, Any], state_file: str) -> str:
    """Format a human-readable summary of the pipeline state."""
    lines = ["Pipeline State Summary:"]
    lines.append(f"  State file: {state_file}")
    lines.append(f"  Pipeline version: {state.get('pipeline_version', 'unknown')}")

    if state.get("start_time"):
        start_time = datetime.fromtimestamp(state["start_time"])
        lines.append(f"  Started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}")

    lines.append("\nSteps:")
    for name, info in state["steps"].items():
        status_symbol = {
            "completed": "✓",
            "failed": "✗",
            "running": "→",
            "finalizing": "◐",
            "pending": "·",
        }.get(info.status, "?")

        duration_str = f" ({info.duration:.1f}s)" if info.duration else ""
        lines.append(f"  {status_symbol} {name}{duration_str}")

        if info.status == "failed" and info.error:
            lines.append(f"    Error: {info.error}")

    return "\n".join(lines)

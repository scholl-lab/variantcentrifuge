"""Helpers for building safe run metadata export views."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any

from ..pipeline_core import PipelineContext

SAFE_CONFIG_KEYS = {
    "reference",
    "genome_build",
    "pipeline_version",
    "threads",
    "inheritance_mode",
    "use_chunked_processing",
    "chunk_size",
    "max_workers",
    "scoring_config_path",
    "preset",
    "filter",
    "filters",
    "final_filter",
    "gzip_intermediates",
    "stats_config",
}


def _sanitize_text_identifiers(
    text: str, mapping: dict[str, str], redact_unmapped: bool = False
) -> str:
    """Replace original identifiers with pseudonyms or redactions in text strings."""
    if not text:
        return text
    if not mapping:
        return "<redacted>" if redact_unmapped else text
    ids = [k for k in sorted(mapping.keys(), key=len, reverse=True) if k]
    if not ids:
        return text
    pattern = re.compile(
        r"(?<![A-Za-z0-9])(" + "|".join(re.escape(i) for i in ids) + r")(?![A-Za-z0-9])"
    )
    return pattern.sub(lambda m: mapping[m.group(1)], text)


def collect_metadata_rows(
    context: PipelineContext,
    get_version_fn: Callable[[str], str],
) -> list[tuple[str, Any]]:
    """Collect parameter/value rows for run metadata TSV and Excel reports.

    When pseudonymization is enabled, sensitive sample identifiers and lookup mappings
    are stripped or replaced to avoid re-identification.
    """
    cfg = context.config
    is_pseudonymized = bool(cfg.get("pseudonymize"))
    mapping: dict[str, str] = {}
    if is_pseudonymized:
        raw_map = cfg.get("pseudonymization_mapping")
        if isinstance(raw_map, dict):
            mapping = {str(k): str(v) for k, v in raw_map.items()}

    raw_preset = cfg.get("preset_name") or cfg.get("preset", "")
    if isinstance(raw_preset, (list, tuple, set)):
        preset_str = ", ".join(str(p) for p in raw_preset)
    else:
        preset_str = str(raw_preset or "")
    if preset_str and cfg.get("preset_description"):
        preset_str = f"{preset_str} ({cfg.get('preset_description')})"

    samples = context.vcf_samples or []

    def _sample_lbl(name_key: str, idx_key: str) -> str:
        if cfg.get(name_key):
            val = str(cfg[name_key])
            return mapping.get(val, val) if is_pseudonymized else val
        idx = cfg.get(idx_key)
        if idx not in (None, "") and str(idx).isdigit() and int(str(idx)) < len(samples):
            val = samples[int(str(idx))]
            return mapping.get(val, val) if is_pseudonymized else val
        return ""

    t_lbl, n_lbl = (
        _sample_lbl("tumor_sample", "tumor_sample_index"),
        _sample_lbl("normal_sample", "normal_sample_index"),
    )
    tn_str = f"tumor={t_lbl or 'N/A'}, normal={n_lbl or 'N/A'}" if (t_lbl or n_lbl) else ""

    tx_file = cfg.get("transcript_file", "")
    tx_md5 = (
        f"{tx_file} ({hashlib.md5(Path(tx_file).read_bytes()).hexdigest()})"
        if tx_file and Path(tx_file).is_file()
        else ""
    )

    vcf_val = str(cfg.get("vcf_file", ""))
    pheno_val = str(cfg.get("phenotype_file", ""))
    ped_val = str(cfg.get("ped_file", ""))
    gene_file_val = str(cfg.get("gene_file", ""))
    filter_val = str(cfg.get("filters") or cfg.get("filter", ""))
    final_filter_val = str(cfg.get("final_filter", "") or "")

    if is_pseudonymized:
        redact = not mapping
        vcf_val = _sanitize_text_identifiers(
            Path(vcf_val).name if vcf_val else "", mapping, redact_unmapped=redact
        )
        pheno_val = _sanitize_text_identifiers(
            Path(pheno_val).name if pheno_val else "", mapping, redact_unmapped=redact
        )
        ped_val = _sanitize_text_identifiers(
            Path(ped_val).name if ped_val else "", mapping, redact_unmapped=redact
        )
        gene_file_val = _sanitize_text_identifiers(
            Path(gene_file_val).name if gene_file_val else "", mapping, redact_unmapped=redact
        )
        filter_val = _sanitize_text_identifiers(filter_val, mapping, redact_unmapped=redact)
        final_filter_val = _sanitize_text_identifiers(
            final_filter_val, mapping, redact_unmapped=redact
        )
        tn_str = _sanitize_text_identifiers(tn_str, mapping, redact_unmapped=redact)
        if redact and tx_md5:
            tx_md5 = "<redacted>"
        elif tx_md5:
            tx_md5 = _sanitize_text_identifiers(tx_md5, mapping, redact_unmapped=redact)

    ref_val = str(cfg.get("reference", ""))
    scoring_val = str(cfg.get("scoring_config_path", ""))
    if is_pseudonymized:
        ref_val = _sanitize_text_identifiers(ref_val, mapping, redact_unmapped=not mapping)
        scoring_val = _sanitize_text_identifiers(
            Path(scoring_val).name if scoring_val else "", mapping, redact_unmapped=not mapping
        )

    rows: list[tuple[str, Any]] = [
        ("Tool", "variantcentrifuge"),
        ("Version", cfg.get("pipeline_version", "0.5.0")),
        ("Git SHA", cfg.get("git_sha") or get_version_fn("git")),
        ("Run_start_time", context.start_time.isoformat()),
        ("Run_duration_seconds", f"{context.get_execution_time():.1f}"),
        ("VCF_file", vcf_val),
        ("Gene_name", cfg.get("normalized_genes", "")),
        ("Gene_file", gene_file_val),
        ("Phenotype_file", pheno_val),
        ("Pedigree_file", ped_val),
        ("Transcript File (MD5)", tx_md5),
        ("Reference", ref_val),
        ("Filter", filter_val),
        ("Preset", preset_str),
        ("Final Filter", final_filter_val),
        ("Tumor/Normal Samples", tn_str),
        ("Threads", cfg.get("threads", "")),
        ("Inheritance_mode", cfg.get("inheritance_mode", "")),
        ("Scoring_config", scoring_val),
    ]

    if is_pseudonymized:
        identifying_config = {
            "reference",
            "scoring_config_path",
            "filter",
            "filters",
            "final_filter",
            "stats_config",
            "preset",
        }
        # Include only safe allowlisted config entries, sanitized
        for k in sorted(SAFE_CONFIG_KEYS):
            if k in cfg and cfg[k] is not None:
                val = str(cfg[k])
                if k in identifying_config:
                    val = _sanitize_text_identifiers(val, mapping, redact_unmapped=not mapping)
                rows.append((f"config.{k}", val))
    else:
        sensitive_keys = {"pseudonymization_mapping", "pseudonymizer"}
        rows.extend(
            (f"config.{k}", str(v))
            for k, v in cfg.items()
            if v is not None and k not in sensitive_keys
        )

    for tool_name in ["bcftools", "snpEff", "SnpSift", "bedtools"]:
        rows.append((f"tool.{tool_name.lower()}_version", get_version_fn(tool_name)))

    return rows

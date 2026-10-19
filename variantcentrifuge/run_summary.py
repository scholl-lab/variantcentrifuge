"""Consolidated run summary and provenance manifest generation."""

from __future__ import annotations

import datetime
import hashlib
import json
import logging
import os
import platform
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .checkpoint_models import AtomicFileOperation
from .version import __version__

if TYPE_CHECKING:
    from .pipeline_core.context import PipelineContext

logger = logging.getLogger("variantcentrifuge")

SUMMARY_SCHEMA_VERSION = "1.0.0"


def compute_stored_sha256(filepath: str | Path, chunk_size: int = 65536) -> str:
    """Compute streamed SHA-256 checksum of raw stored file bytes on disk.

    Reads exact stored bytes without decompression (e.g. for .gz or .xlsx),
    matching standard command-line sha256sum output.
    """
    path = Path(filepath)
    if not path.is_file():
        raise FileNotFoundError(f"File not found for checksum: {filepath}")

    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def _sanitize_for_json(obj: Any) -> Any:
    """Recursively convert objects to JSON-serializable types."""
    if obj is None or isinstance(obj, (int, float, str, bool)):
        return obj
    if isinstance(obj, (Path,)):
        return str(obj)
    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()
    if isinstance(obj, (list, tuple, set)):
        return [_sanitize_for_json(item) for item in obj]
    if isinstance(obj, dict):
        return {str(k): _sanitize_for_json(v) for k, v in obj.items()}
    return str(obj)


def get_environment_versions() -> dict[str, str]:
    """Collect tool, runtime, and dependency versions."""
    from .utils import get_tool_version

    versions: dict[str, str] = {
        "variantcentrifuge": __version__,
        "python": platform.python_version(),
    }
    for tool in ("bcftools", "snpEff", "SnpSift", "bedtools", "git"):
        ver = get_tool_version(tool)
        if ver and ver != "unknown":
            versions[tool] = ver

    return versions


def collect_tier_counts(context: PipelineContext) -> dict[str, int]:
    """Collect variant and entity counts across filtering tiers."""
    counts: dict[str, int] = {}

    if context.vcf_samples:
        counts["samples"] = len(context.vcf_samples)

    if context.current_dataframe is not None:
        counts["final_variants"] = len(context.current_dataframe)
        if "GENE" in context.current_dataframe.columns:
            genes = context.current_dataframe["GENE"].dropna().unique()
            counts["unique_genes"] = len(genes)
    elif context.final_output_path and context.final_output_path.is_file():
        try:
            with open(context.final_output_path, encoding="utf-8") as f:
                lines = sum(1 for _ in f)
            counts["final_variants"] = max(0, lines - 1)
        except Exception:
            pass

    return counts


def collect_file_manifest(context: PipelineContext) -> dict[str, list[dict[str, Any]]]:
    """Collect input and output file paths, byte sizes, and stored SHA-256 checksums."""
    manifest: dict[str, list[dict[str, Any]]] = {"inputs": [], "outputs": []}

    input_keys = [
        "vcf_file",
        "ped_file",
        "phenotype_file",
        "gene_file",
        "regions_bed",
        "transcript_file",
    ]
    seen_inputs: set[str] = set()
    for key in input_keys:
        p_str = context.config.get(key)
        if p_str and os.path.isfile(p_str) and p_str not in seen_inputs:
            seen_inputs.add(p_str)
            manifest["inputs"].append(
                {
                    "type": key,
                    "path": os.path.abspath(p_str),
                    "size_bytes": os.path.getsize(p_str),
                    "sha256": compute_stored_sha256(p_str),
                }
            )

    seen_outputs: set[str] = set()
    output_files: list[tuple[str, Path | str]] = []
    if context.final_output_path:
        output_files.append(("final_output_tsv", context.final_output_path))
    for rtype, rpath in context.report_paths.items():
        output_files.append((f"report_{rtype}", rpath))
    for key in ("gene_burden_output", "association_output"):
        out = context.config.get(key)
        if out:
            output_files.append((key, out))

    for otype, opath in output_files:
        p_str = str(opath)
        if os.path.isfile(p_str) and p_str not in seen_outputs:
            seen_outputs.add(p_str)
            manifest["outputs"].append(
                {
                    "type": otype,
                    "path": os.path.abspath(p_str),
                    "size_bytes": os.path.getsize(p_str),
                    "sha256": compute_stored_sha256(p_str),
                }
            )

    return manifest


def build_run_summary(
    context: PipelineContext,
    status: str = "success",
    error_message: str | None = None,
) -> dict[str, Any]:
    """Build the comprehensive run_summary manifest dictionary."""
    start_iso = context.start_time.isoformat()
    now_iso = datetime.datetime.now().isoformat()
    elapsed = round(context.get_execution_time(), 2)

    stages_info: list[dict[str, Any]] = [
        {"name": stage_name, "status": "completed"}
        for stage_name in sorted(context.completed_stages)
    ]

    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "pipeline_version": __version__,
        "generated_at": now_iso,
        "run_start_time": start_iso,
        "elapsed_seconds": elapsed,
        "execution_status": status,
        "error_message": error_message,
        "configuration": _sanitize_for_json(context.config),
        "versions": get_environment_versions(),
        "stages": stages_info,
        "tier_counts": collect_tier_counts(context),
        "files": collect_file_manifest(context),
    }


def write_run_summary(
    context: PipelineContext,
    output_path: Path | str | None = None,
    status: str = "success",
    error_message: str | None = None,
) -> Path:
    """Atomically write the consolidated run_summary.json file."""
    if output_path is None:
        target_path = context.workspace.output_dir / "run_summary.json"
    else:
        target_path = Path(output_path)

    summary = build_run_summary(context, status=status, error_message=error_message)

    with (
        AtomicFileOperation(target_path) as temp_path,
        open(temp_path, "w", encoding="utf-8") as f,
    ):
        json.dump(summary, f, indent=2)

    logger.info(f"Published run summary: {target_path}")
    return target_path


def _collect_metadata_rows(context: PipelineContext, get_version_fn: Any) -> list[tuple[str, Any]]:
    """Build (Parameter, Value) tuples for MetadataGenerationStage."""
    cfg = context.config
    raw_preset = cfg.get("preset", "")
    preset_str = (
        ", ".join(str(p) for p in raw_preset)
        if isinstance(raw_preset, list | tuple | set)
        else str(raw_preset or "")
    )
    samples = context.vcf_samples or []

    def _sample_lbl(name_key: str, idx_key: str) -> str:
        if cfg.get(name_key):
            return str(cfg[name_key])
        idx = cfg.get(idx_key)
        return (
            samples[int(str(idx))]
            if idx not in (None, "") and str(idx).isdigit() and int(str(idx)) < len(samples)
            else ""
        )

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
    rows: list[tuple[str, Any]] = [
        ("Tool", "variantcentrifuge"),
        ("Version", cfg.get("pipeline_version", __version__)),
        ("Git SHA", cfg.get("git_sha") or get_version_fn("git")),
        ("Run_start_time", context.start_time.isoformat()),
        ("Run_duration_seconds", f"{context.get_execution_time():.1f}"),
        ("VCF_file", cfg.get("vcf_file", "")),
        ("Gene_name", cfg.get("normalized_genes", "")),
        ("Gene_file", cfg.get("gene_file", "")),
        ("Phenotype_file", cfg.get("phenotype_file", "")),
        ("Pedigree_file", cfg.get("ped_file", "")),
        ("Transcript File (MD5)", tx_md5),
        ("Reference", cfg.get("reference", "")),
        ("Filter", cfg.get("filters") or cfg.get("filter", "")),
        ("Preset", preset_str),
        ("Final Filter", cfg.get("final_filter", "") or ""),
        ("Tumor/Normal Samples", tn_str),
        ("Threads", cfg.get("threads", "")),
        ("Inheritance_mode", cfg.get("inheritance_mode", "")),
        ("Scoring_config", cfg.get("scoring_config_path", "")),
    ]
    rows.extend((f"config.{k}", str(v)) for k, v in cfg.items() if v is not None)
    rows.extend(
        [
            ("tool.bcftools_version", get_version_fn("bcftools")),
            ("tool.snpeff_version", get_version_fn("snpEff")),
            ("tool.snpsift_version", get_version_fn("SnpSift")),
            ("tool.bedtools_version", get_version_fn("bedtools")),
        ]
    )
    return rows

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


def _count_tsv_rows(path: Path) -> int:
    """Count data rows in a TSV or compressed TSV file, skipping header."""
    try:
        if str(path).endswith(".gz"):
            import gzip

            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
                lines = sum(1 for _ in f)
        else:
            with open(path, encoding="utf-8", errors="replace") as f:
                lines = sum(1 for _ in f)
        return max(0, lines - 1)
    except Exception:
        return 0


def collect_tier_counts(context: PipelineContext) -> dict[str, Any]:
    """Collect variant and entity counts across filtering tiers."""
    counts: dict[str, Any] = {}

    if context.vcf_samples:
        counts["samples"] = len(context.vcf_samples)

    # Extracted variants
    extracted_tsv = getattr(context, "extracted_tsv", None)
    if extracted_tsv and Path(extracted_tsv).is_file():
        counts["extracted_variants"] = _count_tsv_rows(Path(extracted_tsv))

    # Final variants and gene counts
    df = context.current_dataframe
    if df is not None:
        counts["final_variants"] = len(df)
        gene_col = next(
            (c for c in ("GENE", "ANN[*].GENE", "ANN_0__GENE") if c in df.columns),
            None,
        )
        if gene_col:
            genes = df[gene_col].dropna().unique()
            counts["final_genes"] = len(genes)
            counts["unique_genes"] = len(genes)
        else:
            counts["final_genes"] = 0
            counts["unique_genes"] = 0
    elif context.final_output_path and context.final_output_path.is_file():
        counts["final_variants"] = _count_tsv_rows(context.final_output_path)
        counts["final_genes"] = 0
        counts["unique_genes"] = 0
    else:
        counts["final_variants"] = 0
        counts["final_genes"] = 0
        counts["unique_genes"] = 0

    if "extracted_variants" not in counts:
        counts["extracted_variants"] = counts["final_variants"]

    # Pre-final-filter variants
    pre_final = getattr(context, "pre_final_filter_variants", None)
    if pre_final is None and hasattr(context, "config"):
        pre_final = context.config.get("pre_final_filter_variants")
    if pre_final is not None:
        counts["pre_final_filter_variants"] = int(pre_final)
    else:
        counts["pre_final_filter_variants"] = counts["final_variants"]

    # Impact tiers
    impact_tiers: dict[str, int] = {}
    if df is not None:
        impact_col = next(
            (
                c
                for c in ("IMPACT", "ANN[0].IMPACT", "ANN_0__IMPACT", "ANN[*].IMPACT")
                if c in df.columns
            ),
            None,
        )
        if impact_col:
            val_counts = df[impact_col].value_counts().to_dict()
            for tier in ("HIGH", "MODERATE", "LOW", "MODIFIER"):
                impact_tiers[tier] = int(val_counts.get(tier, 0))
    counts["impact_tiers"] = impact_tiers

    return counts


def collect_file_manifest(context: PipelineContext) -> dict[str, list[dict[str, Any]]]:
    """Collect input and output file paths, byte sizes, and stored SHA-256 checksums."""
    manifest: dict[str, list[dict[str, Any]]] = {"inputs": [], "outputs": []}
    is_pseudonymized = bool(context.config.get("pseudonymize"))
    mapping: dict[str, str] = {}
    if is_pseudonymized:
        raw_map = context.config.get("pseudonymization_mapping")
        if isinstance(raw_map, dict):
            mapping = {str(k): str(v) for k, v in raw_map.items()}

    def _manifest_path(p_str: str) -> str:
        if not is_pseudonymized:
            return os.path.abspath(p_str)
        from .stages.metadata_helpers import _sanitize_text_identifiers

        return _sanitize_text_identifiers(Path(p_str).name, mapping, redact_unmapped=not mapping)

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
                    "path": _manifest_path(p_str),
                    "size_bytes": os.path.getsize(p_str),
                    "sha256": compute_stored_sha256(p_str),
                }
            )

    seen_outputs: set[str] = set()
    output_files: list[tuple[str, Path | str]] = []
    if context.final_output_path:
        output_files.append(("final_output_tsv", context.final_output_path))
    for rtype, rpath in context.report_paths.items():
        if rtype != "run_summary" and not str(rpath).endswith("run_summary.json"):
            output_files.append((f"report_{rtype}", rpath))
    for key in ("gene_burden_output", "association_output"):
        out = context.config.get(key)
        if out:
            output_files.append((key, out))

    if context.stats_output_tsv:
        output_files.append(("stats_output_file", context.stats_output_tsv))
    elif context.config.get("stats_output_file"):
        stats_out = Path(context.config["stats_output_file"])
        if not stats_out.is_absolute() and context.workspace:
            stats_out = context.workspace.output_dir / stats_out
        output_files.append(("stats_output_file", stats_out))

    if context.workspace and (context.workspace.output_dir / "report" / "summary.json").is_file():
        output_files.append(
            ("report_summary_json", context.workspace.output_dir / "report" / "summary.json")
        )

    for otype, opath in output_files:
        p_str = str(opath)
        if os.path.isfile(p_str) and p_str not in seen_outputs:
            seen_outputs.add(p_str)
            manifest["outputs"].append(
                {
                    "type": otype,
                    "path": _manifest_path(p_str),
                    "size_bytes": os.path.getsize(p_str),
                    "sha256": compute_stored_sha256(p_str),
                }
            )

    return manifest


def _build_safe_configuration(context: PipelineContext) -> dict[str, Any]:
    """Build sanitized configuration dictionary for run_summary.json."""
    if context.config.get("pseudonymize"):
        rows = _collect_metadata_rows(context, lambda _tool: "")
        return {k[7:]: v for k, v in rows if k.startswith("config.")}
    return {
        str(k): _sanitize_for_json(v)
        for k, v in context.config.items()
        if k not in ("pseudonymization_mapping", "pseudonymizer")
    }


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
        {
            "name": stage_name,
            "status": "completed",
            "duration_seconds": context.stage_durations.get(stage_name),
        }
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
        "configuration": _build_safe_configuration(context),
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


from .stages.metadata_helpers import collect_metadata_rows as _collect_metadata_rows

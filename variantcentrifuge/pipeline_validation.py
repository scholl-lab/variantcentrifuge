"""Output validation helpers for the pipeline runner."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .pipeline_core import PipelineContext


def validate_requested_outputs(context: PipelineContext) -> None:
    """Ensure requested user-facing outputs exist before reporting success."""
    output_file = context.config.get("output_file")
    if output_file not in (None, "stdout", "-") and (
        not context.final_output_path or not context.final_output_path.exists()
    ):
        raise RuntimeError("Requested TSV output was not written")

    if context.config.get("perform_association"):
        assoc = context.config.get("association_output")
        if not assoc or not Path(assoc).exists():
            raise RuntimeError("Requested association output was not written")

    if context.config.get("perform_gene_burden"):
        burden = context.config.get("gene_burden_output")
        if not burden or not Path(burden).exists():
            raise RuntimeError("Requested gene burden output was not written")

    if context.config.get("xlsx") or context.config.get("excel"):
        excel = context.report_paths.get("excel")
        if not excel or not Path(excel).exists():
            raise RuntimeError("Requested Excel output was not written")

    if context.config.get("html_report") or context.config.get("html"):
        html = context.report_paths.get("html")
        if not html or not Path(html).exists():
            raise RuntimeError("Requested HTML report output was not written")


def resolve_threads(raw_threads: object, config: dict | None = None) -> int:
    """Resolve thread setting to an integer count."""
    if str(raw_threads).lower() == "auto":
        from .memory.resource_manager import ResourceManager

        rm = ResourceManager(config=config or {})
        return int(rm.cpu_cores)
    try:
        return int(str(raw_threads))
    except (ValueError, TypeError):
        return 1

"""Chunk TSV merging and cleanup for ParallelCompleteProcessingStage (split for size)."""

from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path

from ..pipeline_core import PipelineContext

logger = logging.getLogger(__name__)


class ChunkTsvMergeMixin:
    """Merge per-chunk extracted TSVs into one file and remove chunk files."""

    def _merge_tsv_outputs(self, context: PipelineContext, chunk_tsvs: list[Path]) -> Path:
        """Merge TSV chunks into a single file with optimized compression handling."""
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.tsv"
        )
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression:
            output_tsv = Path(str(output_tsv) + ".gz")

        expected_header = self._build_expected_header(context)

        if not chunk_tsvs:
            hdr = (expected_header or "") + "\n"
            if use_compression:
                with gzip.open(output_tsv, "wt", compresslevel=1) as f:
                    f.write(hdr)
            else:
                output_tsv.write_text(hdr, encoding="utf-8")
            logger.info(f"Created schema-valid header-only TSV for empty chunks: {output_tsv}")
            return output_tsv

        if len(chunk_tsvs) == 1:
            if expected_header is not None:
                open_fn = gzip.open if str(chunk_tsvs[0]).endswith(".gz") else open
                with open_fn(str(chunk_tsvs[0]), "rt") as in_fh:
                    hdr = in_fh.readline().strip()
                if hdr != expected_header.strip():
                    raise ValueError(
                        f"Incompatible chunk header in {chunk_tsvs[0]}: "
                        f"expected '{expected_header.strip()}', got '{hdr}'"
                    )
            shutil.move(str(chunk_tsvs[0]), str(output_tsv))
            logger.info(f"Moved single TSV chunk to {output_tsv}")
        else:
            self._merge_tsvs_traditional(chunk_tsvs, output_tsv, expected_header=expected_header)
            logger.info(f"Merged {len(chunk_tsvs)} TSV chunks into {output_tsv}")

        return output_tsv

    def _merge_tsvs_traditional(
        self,
        chunk_tsvs: list[Path],
        output_tsv: Path,
        expected_header: str | None = None,
    ) -> None:
        with (
            gzip.open(output_tsv, "wt", compresslevel=1)
            if str(output_tsv).endswith(".gz")
            else open(output_tsv, "w")
        ) as out_fh:
            first_header: str | None = None
            for chunk_tsv in sorted(chunk_tsvs):
                with (
                    gzip.open(chunk_tsv, "rt")
                    if str(chunk_tsv).endswith(".gz")
                    else open(chunk_tsv)
                ) as in_fh:
                    header = in_fh.readline()
                    hdr_clean = header.strip()
                    if expected_header is not None and hdr_clean != expected_header.strip():
                        raise ValueError(
                            f"Incompatible chunk header in {chunk_tsv}: "
                            f"expected '{expected_header.strip()}', got '{hdr_clean}'"
                        )
                    if first_header is None:
                        first_header = hdr_clean
                        out_fh.write(header)
                    elif hdr_clean != first_header:
                        raise ValueError(
                            f"Incompatible chunk header in {chunk_tsv}: "
                            f"expected '{first_header}', got '{hdr_clean}'"
                        )
                    shutil.copyfileobj(in_fh, out_fh)

    def _cleanup_chunks(
        self, bed_chunks: list[Path], tsv_chunks: list[Path], context: PipelineContext
    ) -> None:
        """Clean up temporary chunk files."""
        if not context.config.get("keep_intermediates", False):
            for chunk in bed_chunks + tsv_chunks:
                if chunk.exists():
                    chunk.unlink()
                meta = Path(f"{chunk}.meta")
                if meta.exists():
                    meta.unlink()

            if bed_chunks and bed_chunks[0].parent.name == "chunks":
                chunk_dir = bed_chunks[0].parent
                if chunk_dir.exists() and not any(chunk_dir.iterdir()):
                    chunk_dir.rmdir()

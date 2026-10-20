"""Chunked variant analysis stage."""

from __future__ import annotations

import csv
import logging
from collections.abc import Iterator
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Any, cast

import pandas as pd

from ..dataframe_optimizer import (
    load_optimized_dataframe,
    rename_invalid_identifiers,
    should_use_memory_passthrough,
)
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from .chunk_helpers import (
    TSV_READ_KWARGS,
    apply_chunk_annotations,
    apply_chunk_inheritance_analysis,
    arrow_tsv_batches,
    ensure_sorted_by_gene,
    read_tsv_in_gene_chunks,
    serial_tsv_header,
)
from .genotype_matrix_builder import _STANDARD_COLUMN_ALIASES, _create_standard_column_aliases

logger = logging.getLogger(__name__)

_ROW_ORDER_COL = "_vc_input_row"
_OTHER_GENE_COLUMNS = ("Gene", "gene", "SYMBOL")


def _find_gene_column(header: list[str]) -> str | None:
    """Return the raw header column that the non-chunked loader would expose as ``GENE``.

    Same resolution as DataFrameLoadingStage: an explicit GENE column, else the first alias
    source in ``_create_standard_column_aliases`` order; other gene-like names last.
    """
    _, sanitized = rename_invalid_identifiers(pd.DataFrame(columns=header))
    raw_by_name = {sanitized.get(col, col): col for col in header}
    aliases = [src for src, short in _STANDARD_COLUMN_ALIASES.items() if short == "GENE"]
    for name in ("GENE", *aliases, *_OTHER_GENE_COLUMNS):
        if name in raw_by_name:
            return raw_by_name[name]
    return None


def _write_row_order_copy(input_file: Path, output_file: Path, chunksize: int) -> None:
    """Write the rows the serial loader parses, prefixed by their parsed position.

    Mirrors load_optimized_dataframe's two tiers in bounded batches: PyArrow first; on a
    short row or duplicate header, the pandas C engine with the same settings. Blank lines,
    a BOM, gzip and malformed rows are therefore handled exactly as in the serial path.
    """
    header = serial_tsv_header(input_file)
    if len(set(header)) == len(header):
        try:
            _write_indexed(arrow_tsv_batches(input_file, header), header, output_file)
            return
        except ValueError as e:  # pyarrow.ArrowInvalid, e.g. a short row
            logger.debug(f"PyArrow parse failed ({e}); using the pandas C engine like serial")
    # pandas mis-parses malformed rows in tiny chunks; the floor keeps it batch-invariant.
    batches = pd.read_csv(input_file, chunksize=max(chunksize, 1000), **TSV_READ_KWARGS)
    _write_indexed(batches, _read_header(input_file), output_file)


def _write_indexed(batches: Iterator[pd.DataFrame], header: list[str], output_file: Path) -> None:
    with open(output_file, "w") as dst:
        pd.DataFrame(columns=[_ROW_ORDER_COL, *header]).to_csv(dst, sep="\t", index=False)
        start = 0
        for part in batches:
            part.insert(0, _ROW_ORDER_COL, range(start, start + len(part)))
            start += len(part)
            part.to_csv(dst, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)


def _read_header(path: Path) -> list[str]:
    """Column names as the serial reader parses them (BOM stripped, gzip handled)."""
    return list(pd.read_csv(path, nrows=0, **TSV_READ_KWARGS).columns)


def _load_chunk_artifact(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load a chunked-analysis artifact into the serial loader's schema and row order."""
    df, rename_map = load_optimized_dataframe(
        str(path),
        sep="\t",
        compression="gzip" if str(path).endswith(".gz") else None,
        quoting=csv.QUOTE_MINIMAL,
        on_bad_lines="warn",
    )
    _create_standard_column_aliases(df)
    return _restore_input_order(df), rename_map


def _restore_input_order(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows in input order (as the non-chunked path keeps them) without the helper column."""
    if _ROW_ORDER_COL not in df.columns:
        return df
    order = pd.to_numeric(df[_ROW_ORDER_COL]).to_numpy()
    df = df.drop(columns=_ROW_ORDER_COL).iloc[order.argsort(kind="stable")]
    return df.reset_index(drop=True)


def _normalize_chunk(chunk_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    """Give a raw chunk the column names the non-chunked DataFrame loader produces."""
    chunk_df, rename_map = rename_invalid_identifiers(chunk_df)
    _create_standard_column_aliases(chunk_df)
    return chunk_df, rename_map


def _align_chunk_columns(chunk: pd.DataFrame, expected_columns: list[str]) -> pd.DataFrame:
    """Reorder a chunk to the first chunk's columns; refuse to drop or invent columns."""
    if set(chunk.columns) != set(expected_columns):
        extra = sorted(set(chunk.columns) - set(expected_columns))
        missing = sorted(set(expected_columns) - set(chunk.columns))
        raise ValueError(f"Chunk columns differ from first chunk: extra={extra}, missing={missing}")
    return chunk.reindex(columns=expected_columns)


class ChunkedAnalysisStage(Stage):
    """Process large files in chunks with all analysis steps."""

    def __init__(self, max_variants_per_gene: int = 10000):
        super().__init__()
        self.max_variants_per_gene = max_variants_per_gene

    @property
    def name(self) -> str:
        return "chunked_analysis"

    @property
    def checkpoint_resume_policy(self) -> str:
        return "restore"

    @property
    def description(self) -> str:
        return "Process large files in memory-efficient chunks"

    @property
    def dependencies(self) -> set[str]:
        return set()

    @property
    def soft_dependencies(self) -> set[str]:
        return {
            "field_extraction",
            "phenotype_integration",
            "genotype_replacement",
            "parallel_complete_processing",
            "extra_column_removal",
            "dataframe_loading",
            "data_sorting",
            "annotation_config_loading",
            "scoring_config_loading",
            "pedigree_loading",
        }

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Restore chunked DataFrame analysis state from its durable TSV artifact."""
        if not context.config.get("use_chunked_processing"):
            # The completed run was a no-op (see _process); it created no artifact.
            logger.info("chunked_analysis was a no-op in the checkpointed run; nothing to restore")
            return context

        chunked_tsv = context.chunked_analysis_tsv
        if not chunked_tsv:
            base = context.workspace.get_intermediate_path("chunked_analysis_results.tsv")
            gz = Path(str(base) + ".gz")
            chunked_tsv = gz if gz.exists() else base

        if not chunked_tsv or not Path(chunked_tsv).exists():
            input_file = self._find_input_file(context)
            if _find_gene_column(_read_header(input_file)) is None:
                # The checkpointed run took the non-chunked fallback; replay it.
                return self._load_unchunked(context, input_file)
            raise RuntimeError(
                "Cannot restore chunked_analysis from checkpoint: "
                f"chunked output missing: {chunked_tsv}"
            )

        df, rename_map = _load_chunk_artifact(Path(chunked_tsv))
        context.current_dataframe = df
        context.column_rename_map = rename_map
        context.chunked_analysis_tsv = Path(chunked_tsv)
        context.config["dataframe_chunked_analysis_complete"] = True
        context.config["use_chunked_processing"] = True
        logger.info(
            "Restored %d variants from chunked analysis checkpoint artifact: %s",
            len(df),
            chunked_tsv,
        )
        return context

    def _setup_inheritance_config(self, context: PipelineContext) -> None:
        """Set up inheritance analysis configuration for chunked processing."""
        should_calculate_inheritance = self._should_calculate_inheritance(context)

        if should_calculate_inheritance:
            logger.info("Setting up inheritance analysis configuration for chunked processing")
            inheritance_mode = context.config.get("inheritance_mode", "simple")
            preserve_details = self._check_if_scoring_needs_details(context)
            if preserve_details:
                logger.debug("Preserving Inheritance_Details for chunked processing scoring")

            context.config["inheritance_analysis_config"] = {
                "inheritance_mode": inheritance_mode,
                "vcf_samples": context.vcf_samples,
                "pedigree_data": context.pedigree_data,
                "use_vectorized_comp_het": not context.config.get("no_vectorized_comp_het", False),
                "preserve_details_for_scoring": preserve_details,
            }

            logger.debug(
                f"Inheritance config set up for chunked processing: mode={inheritance_mode}, "
                f"samples={len(context.vcf_samples or [])}, "
                f"pedigree={len(context.pedigree_data or {})}"
            )
        else:
            logger.debug("Inheritance analysis not needed for chunked processing")

    def _should_calculate_inheritance(self, context: PipelineContext) -> bool:
        """Check if inheritance analysis should be calculated."""
        has_ped_file = context.pedigree_data is not None and len(context.pedigree_data) > 0
        has_inheritance_mode = context.config.get("inheritance_mode")
        has_calculate_inheritance_config = context.config.get("calculate_inheritance", False)
        requires_inheritance_for_scoring = self._check_if_scoring_needs_details(context)

        should_calculate = (
            has_ped_file
            or has_inheritance_mode
            or has_calculate_inheritance_config
            or requires_inheritance_for_scoring
        )
        return bool(should_calculate)

    def _check_if_scoring_needs_details(self, context: PipelineContext) -> bool:
        """Check if scoring configuration requires Inheritance_Details column."""
        if not context.scoring_config:
            return False

        if "formulas" in context.scoring_config:
            for formula_dict in context.scoring_config["formulas"]:
                for formula_expr in formula_dict.values():
                    if isinstance(formula_expr, str) and any(
                        var in formula_expr
                        for var in [
                            "Inheritance_Pattern",
                            "Inheritance_Details",
                            "details",
                            "segregation_p_value",
                        ]
                    ):
                        return True

        if "variables" in context.scoring_config:
            variables = context.scoring_config["variables"]
            if any(col in variables for col in ["Inheritance_Pattern", "Inheritance_Details"]):
                return True

        return False

    @staticmethod
    def _find_input_file(context: PipelineContext) -> Path:
        if context.extra_columns_removed_tsv:
            return Path(context.extra_columns_removed_tsv)
        if context.phenotypes_added_tsv:
            return Path(context.phenotypes_added_tsv)
        if context.genotype_replaced_tsv:
            return Path(context.genotype_replaced_tsv)
        if context.extracted_tsv:
            return Path(context.extracted_tsv)
        if context.data:
            return Path(context.data)
        raise RuntimeError("No TSV file available for chunked processing")

    @staticmethod
    def _load_unchunked(context: PipelineContext, input_file: Path) -> PipelineContext:
        """Load the whole TSV exactly like DataFrameLoadingStage and disable chunking."""
        df, rename_map = load_optimized_dataframe(
            str(input_file),
            sep="\t",
            compression="gzip" if str(input_file).endswith(".gz") else None,
            on_bad_lines="warn",
        )
        _create_standard_column_aliases(df)
        context.current_dataframe = df
        context.column_rename_map = rename_map
        if should_use_memory_passthrough(df):
            context.variants_df = df
        context.config["use_chunked_processing"] = False
        logger.info(f"Loaded {len(df)} variants without chunking from {input_file}")
        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in chunks."""
        if not context.config.get("use_chunked_processing"):
            logger.debug("Chunked processing not enabled")
            return context

        input_file = self._find_input_file(context)
        gene_col_name = _find_gene_column(_read_header(input_file))
        if not gene_col_name:
            logger.warning("No GENE column found, falling back to non-chunked processing")
            return self._load_unchunked(context, input_file)

        self._setup_inheritance_config(context)

        inheritance_config = context.config.get("inheritance_analysis_config", {})
        vcf_samples = inheritance_config.get("vcf_samples", [])
        num_samples = len(vcf_samples) if vcf_samples else 1
        total_variants = 100000

        rm = context.resource_manager
        if rm is None:
            from ..memory import ResourceManager

            rm = ResourceManager(config=context.config)
        chunk_size = rm.auto_chunk_size(total_variants, num_samples)
        logger.info(f"Processing {input_file} in chunks of {chunk_size} variants")

        ordered_file = context.workspace.get_intermediate_path("chunked_analysis_input.tsv")
        _write_row_order_copy(input_file, ordered_file, chunk_size)
        # ensure_sorted_by_gene reuses an existing sorted file; never reuse a stale one.
        ordered_file.with_suffix(".sorted.tsv").unlink(missing_ok=True)
        sorted_file = self._ensure_sorted_by_gene(ordered_file, gene_col_name, context)

        threads_config = context.config.get("threads", 1)
        memory_per_chunk_gb = rm.estimate_memory(chunk_size, num_samples)
        memory_capped_workers = rm.auto_workers(
            task_count=max(1, threads_config),
            memory_per_task_gb=memory_per_chunk_gb,
        )
        max_workers = min(threads_config, memory_capped_workers)
        parallel_chunks = context.config.get("parallel_chunks", False) or max_workers > 1

        chunked_output_path = context.workspace.get_intermediate_path(
            "chunked_analysis_results.tsv"
        )
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression:
            chunked_output_path = Path(str(chunked_output_path) + ".gz")
            compression: str | None = "gzip"
        else:
            compression = None

        if parallel_chunks and max_workers > 1:
            logger.info(f"Processing chunks in parallel with {max_workers} workers")
            chunks_written, _ = self._process_chunks_parallel(
                sorted_file,
                gene_col_name,
                chunk_size,
                context,
                max_workers,
                chunked_output_path,
                compression,
            )
        else:
            logger.info("Processing chunks sequentially")
            chunks_written, _ = self._process_chunks_sequential(
                sorted_file,
                gene_col_name,
                chunk_size,
                context,
                chunked_output_path,
                compression,
            )

        if chunks_written == 0:
            logger.warning("No chunks with variants were processed - creating empty DataFrame")
            empty_df = pd.DataFrame(columns=_read_header(sorted_file))
            empty_df.to_csv(
                chunked_output_path, sep="\t", index=False, compression=cast(Any, compression)
            )
        # Empty or not, the artifact is loaded exactly as a checkpoint restore loads it.
        df, rename_map = _load_chunk_artifact(chunked_output_path)
        context.current_dataframe = df
        context.column_rename_map = rename_map

        context.config["dataframe_chunked_analysis_complete"] = True
        context.config.pop("chunked_processing_complete", None)
        context.chunked_analysis_tsv = chunked_output_path
        logger.info(
            f"Chunked processing completed: {len(context.current_dataframe)} total variants"
        )
        logger.info(f"Chunked analysis results written to: {chunked_output_path}")

        return context

    def _process_chunks_sequential(
        self,
        sorted_file: Path,
        gene_col_name: str,
        chunk_size: int,
        context: PipelineContext,
        output_path: Path,
        compression: str | None,
    ) -> tuple[int, int]:
        """Process chunks sequentially and stream directly to disk."""
        chunks_written = 0
        total_variants = 0
        expected_columns: list[str] | None = None
        for chunk_num, chunk_df in enumerate(
            self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size, context)
        ):
            logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")
            processed_chunk = self._process_single_chunk(chunk_df, context, chunk_num)
            if processed_chunk is not None and not processed_chunk.empty:
                if chunks_written == 0:
                    expected_columns = list(processed_chunk.columns)
                    processed_chunk.to_csv(
                        output_path,
                        sep="\t",
                        index=False,
                        mode="w",
                        compression=cast(Any, compression),
                    )
                else:
                    if expected_columns is not None:
                        processed_chunk = _align_chunk_columns(processed_chunk, expected_columns)
                    processed_chunk.to_csv(
                        output_path,
                        sep="\t",
                        index=False,
                        mode="a",
                        header=False,
                        compression=cast(Any, compression),
                    )
                chunks_written += 1
                total_variants += len(processed_chunk)
        return chunks_written, total_variants

    def _process_chunks_parallel(
        self,
        sorted_file: Path,
        gene_col_name: str,
        chunk_size: int,
        context: PipelineContext,
        max_workers: int,
        output_path: Path,
        compression: str | None,
    ) -> tuple[int, int]:
        """Process chunks in parallel with bounded backlog and stream directly to disk."""
        chunk_iter = self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size, context)
        buffered_results: dict[int, pd.DataFrame] = {}

        max_pending = max(max_workers * 2, 4)
        chunks_written = 0
        total_variants = 0
        next_write_idx = 0
        expected_columns: list[str] | None = None

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_index: dict[Any, int] = {}
            chunk_index = 0
            has_more_chunks = True

            while has_more_chunks or future_to_index or buffered_results:
                while (
                    has_more_chunks and (len(future_to_index) + len(buffered_results)) < max_pending
                ):
                    try:
                        chunk_df = next(chunk_iter)
                        future = executor.submit(
                            self._process_single_chunk, chunk_df, context, chunk_index
                        )
                        future_to_index[future] = chunk_index
                        chunk_index += 1
                    except StopIteration:
                        has_more_chunks = False
                        break

                if not future_to_index and not buffered_results:
                    break

                if future_to_index:
                    done, _ = wait(future_to_index.keys(), return_when=FIRST_COMPLETED)
                    for future in done:
                        idx = future_to_index.pop(future)
                        try:
                            processed_chunk = future.result()
                            buffered_results[idx] = processed_chunk
                            logger.info(f"Completed processing chunk {idx + 1}")
                        except Exception as e:
                            logger.error(f"Error processing chunk {idx + 1}: {e}")
                            raise

                while next_write_idx in buffered_results:
                    chunk = buffered_results.pop(next_write_idx)
                    next_write_idx += 1
                    if chunk is not None and not chunk.empty:
                        if chunks_written == 0:
                            expected_columns = list(chunk.columns)
                            chunk.to_csv(
                                output_path,
                                sep="\t",
                                index=False,
                                mode="w",
                                compression=cast(Any, compression),
                            )
                        else:
                            if expected_columns is not None:
                                chunk = _align_chunk_columns(chunk, expected_columns)
                            chunk.to_csv(
                                output_path,
                                sep="\t",
                                index=False,
                                mode="a",
                                header=False,
                                compression=cast(Any, compression),
                            )
                        chunks_written += 1
                        total_variants += len(chunk)

        return chunks_written, total_variants

    def _process_single_chunk(
        self, chunk_df: pd.DataFrame, context: PipelineContext, chunk_num: int
    ) -> pd.DataFrame:
        """Process a single chunk with all analysis steps."""
        logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")

        if chunk_df.empty:
            logger.warning(f"Chunk {chunk_num + 1} is empty, skipping processing")
            return chunk_df

        chunk_df, rename_map = _normalize_chunk(chunk_df)

        if context.annotation_configs:
            chunk_df = self._apply_chunk_annotations(chunk_df, context)

        inheritance_config = context.config.get("inheritance_analysis_config")
        if inheritance_config and not context.config.get("skip_inheritance", False):
            logger.debug(f"Applying inheritance analysis to chunk {chunk_num + 1}")
            chunk_df = self._apply_chunk_inheritance_analysis(chunk_df, context)
            logger.debug(f"Chunk {chunk_num + 1} after inheritance: {list(chunk_df.columns)}")
        else:
            logger.warning(
                f"Inheritance analysis not applied to chunk {chunk_num + 1}: "
                f"config={inheritance_config is not None}, "
                f"skip={context.config.get('skip_inheritance', False)}"
            )

        if context.scoring_config:
            logger.debug(f"Applying scoring to chunk {chunk_num + 1}")
            chunk_df = apply_scoring(chunk_df, context.scoring_config)

        logger.debug(f"Chunk {chunk_num + 1} final size: {len(chunk_df)} variants")
        # Write raw names back; the final reload re-sanitizes them like the non-chunked loader.
        return chunk_df.rename(columns={v: k for k, v in rename_map.items()})

    def _ensure_sorted_by_gene(
        self, input_file: Path, gene_col: str, context: PipelineContext
    ) -> Path:
        return ensure_sorted_by_gene(input_file, gene_col, context)

    def _read_tsv_in_gene_chunks(
        self,
        filepath: Path,
        gene_column: str,
        chunksize: int,
        context: PipelineContext | None = None,
    ) -> Iterator[pd.DataFrame]:
        return read_tsv_in_gene_chunks(
            filepath, gene_column, chunksize, context, self.max_variants_per_gene
        )

    def _apply_chunk_annotations(
        self, chunk_df: pd.DataFrame, context: PipelineContext
    ) -> pd.DataFrame:
        return apply_chunk_annotations(chunk_df, context)

    def _apply_chunk_inheritance_analysis(
        self, chunk_df: pd.DataFrame, context: PipelineContext
    ) -> pd.DataFrame:
        return apply_chunk_inheritance_analysis(chunk_df, context)

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        for attr in [
            "extra_columns_removed_tsv",
            "phenotypes_added_tsv",
            "genotype_replaced_tsv",
            "extracted_tsv",
        ]:
            if hasattr(context, attr):
                tsv_file = getattr(context, attr)
                if tsv_file and tsv_file.exists():
                    return [tsv_file]
        if hasattr(context, "data") and context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        # Only a real chunked run creates the artifact; a no-op/fallback records none.
        return [Path(context.chunked_analysis_tsv)] if context.chunked_analysis_tsv else []

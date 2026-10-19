"""Chunked variant analysis stage."""

from __future__ import annotations

import csv
import gzip
import logging
from collections.abc import Iterator
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Any, cast

import pandas as pd

from ..dataframe_optimizer import load_optimized_dataframe
from ..pipeline_core import PipelineContext, Stage
from ..scoring import apply_scoring
from .chunk_helpers import (
    apply_chunk_annotations,
    apply_chunk_inheritance_analysis,
    calculate_inheritance_safe_chunk_size,
    ensure_sorted_by_gene,
    read_tsv_in_gene_chunks,
)

logger = logging.getLogger(__name__)


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
        chunked_tsv = context.chunked_analysis_tsv
        if not chunked_tsv:
            base = context.workspace.get_intermediate_path("chunked_analysis_results.tsv")
            gz = Path(str(base) + ".gz")
            chunked_tsv = gz if gz.exists() else base

        if not chunked_tsv or not Path(chunked_tsv).exists():
            raise RuntimeError(
                "Cannot restore chunked_analysis from checkpoint: "
                f"chunked output missing: {chunked_tsv}"
            )

        compression = "gzip" if str(chunked_tsv).endswith(".gz") else None
        df, rename_map = load_optimized_dataframe(
            str(chunked_tsv),
            sep="\t",
            compression=compression,
            quoting=csv.QUOTE_MINIMAL,
            on_bad_lines="warn",
        )
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

    def _has_configs(self) -> bool:
        return True

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

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in chunks."""
        if not context.config.get("use_chunked_processing"):
            logger.debug("Chunked processing not enabled")
            return context

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

        if context.extra_columns_removed_tsv:
            input_file = context.extra_columns_removed_tsv
        elif context.phenotypes_added_tsv:
            input_file = context.phenotypes_added_tsv
        elif context.genotype_replaced_tsv:
            input_file = context.genotype_replaced_tsv
        elif context.extracted_tsv:
            input_file = context.extracted_tsv
        elif context.data:
            input_file = Path(context.data)
        else:
            raise RuntimeError("No TSV file available for chunked processing")

        logger.info(f"Processing {input_file} in chunks of {chunk_size} variants")

        gene_col_name = None
        open_fn: Any = gzip.open if str(input_file).endswith(".gz") else open
        with open_fn(input_file, "rt") as f:
            header = f.readline().strip().split("\t")
            for gene_col in ["GENE", "Gene", "gene", "SYMBOL", "ANN[*].GENE", "EFF[*].GENE"]:
                if gene_col in header:
                    gene_col_name = gene_col
                    break

        if not gene_col_name:
            logger.warning("No GENE column found, falling back to non-chunked processing")
            return context

        sorted_file = self._ensure_sorted_by_gene(input_file, gene_col_name, context)

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

        open_func = gzip.open if str(sorted_file).endswith(".gz") else open
        with open_func(sorted_file, "rt") as f:
            first_line = f.readline().strip()
            raw_header = first_line.split("\t") if first_line else []

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
            empty_df = pd.DataFrame(columns=raw_header)
            empty_df.to_csv(
                chunked_output_path, sep="\t", index=False, compression=cast(Any, compression)
            )
            context.current_dataframe = empty_df
        else:
            df, rename_map = load_optimized_dataframe(
                str(chunked_output_path),
                sep="\t",
                compression=compression,
                quoting=csv.QUOTE_MINIMAL,
                on_bad_lines="warn",
            )
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
        for chunk_num, chunk_df in enumerate(
            self._read_tsv_in_gene_chunks(sorted_file, gene_col_name, chunk_size, context)
        ):
            logger.debug(f"Processing chunk {chunk_num + 1} with {len(chunk_df)} variants")
            processed_chunk = self._process_single_chunk(chunk_df, context, chunk_num)
            if processed_chunk is not None and not processed_chunk.empty:
                if chunks_written == 0:
                    processed_chunk.to_csv(
                        output_path,
                        sep="\t",
                        index=False,
                        mode="w",
                        compression=cast(Any, compression),
                    )
                else:
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
                            chunk.to_csv(
                                output_path,
                                sep="\t",
                                index=False,
                                mode="w",
                                compression=cast(Any, compression),
                            )
                        else:
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
        return chunk_df

    def _ensure_sorted_by_gene(
        self, input_file: Path, gene_col: str, context: PipelineContext
    ) -> Path:
        return ensure_sorted_by_gene(input_file, gene_col, context)

    def _calculate_inheritance_safe_chunk_size(
        self, sample_count: int, context: PipelineContext | None = None
    ) -> int:
        return calculate_inheritance_safe_chunk_size(sample_count, context)

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
        return []

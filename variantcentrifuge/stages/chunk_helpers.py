"""Helper functions for chunked variant analysis."""

from __future__ import annotations

import contextlib
import csv
import gzip
import logging
import os
import shutil
import subprocess
import tempfile
import time
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd

from ..gt_reconstruction import _copy_gt_columns_to_samples, _find_per_sample_gt_columns
from ..pipeline_core import PipelineContext
from .sample_column_helpers import (
    cleanup_sample_columns,
    create_sample_columns_from_gt_intelligent,
)

logger = logging.getLogger(__name__)

# Values parse as in the serial loader: only "" is missing (REF=NA survives), quotes are literal.
TSV_READ_KWARGS: dict[str, Any] = {
    "sep": "\t",
    "dtype": str,
    "keep_default_na": False,
    "na_values": [""],
    "quoting": csv.QUOTE_NONE,
    "on_bad_lines": "warn",
    "low_memory": False,
}


_HINT = "Fix the input file or rerun with --no-chunked-processing."


def _tsv_lines(path: Path) -> Iterator[tuple[int, str]]:
    """Yield (line number, line) for non-blank lines, as the serial loader sees them.

    Newlines are universal (CRLF/CR), a UTF-8 BOM on line 1 is dropped, and lines that are
    empty or contain only ASCII spaces are skipped. Any other line is a record, so a line of
    tabs, form feeds or Unicode spaces reaches width validation. Parity exception: in a
    one-column file a spaces-only record is kept by the serial loader but skipped here.
    """
    opener: Any = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if lineno == 1:
                line = line.removeprefix("\ufeff")
            if line.strip(" "):
                yield lineno, line


def _check_no_nul(path: Path, lineno: int, line: str) -> None:
    if "\x00" in line:
        raise ValueError(
            f"Chunked processing: {path} line {lineno} contains a NUL character. {_HINT}"
        )


def read_tsv_header(path: Path) -> list[str]:
    """Header fields (first non-blank line), verbatim apart from BOM/newline.

    Chunked processing requires a header with non-empty field names and no NUL characters.
    """
    lineno, header = next(_tsv_lines(path), (0, ""))
    if not header:
        raise ValueError(f"Chunked processing: {path} has no header line. {_HINT}")
    _check_no_nul(path, lineno, header)
    fields = header.split("\t")
    if "" in fields:
        raise ValueError(
            f"Chunked processing: {path} header (line {lineno}) has an empty column name "
            f"at position {fields.index('') + 1}. {_HINT}"
        )
    return fields


def write_row_indexed_copy(input_file: Path, output_file: Path, index_col: str) -> None:
    """Copy a TSV with each data row prefixed by its 0-based position.

    Fast single text pass (plain or gzip). Chunked processing requires well-formed input:
    a valid header (see read_tsv_header) and data rows with exactly the header's field
    count and no NUL characters. Unlike the serial loader, which skips or pads malformed
    rows, any violation fails here with its line number.
    """
    fields = read_tsv_header(input_file)
    lines = _tsv_lines(input_file)
    next(lines)  # header, validated above
    with open(output_file, "w", encoding="utf-8") as dst:
        dst.write(index_col + "\t" + "\t".join(fields) + "\n")
        for row, (lineno, line) in enumerate(lines):
            if line.count("\t") != len(fields) - 1:
                raise ValueError(
                    f"Chunked processing: {input_file} line {lineno} has "
                    f"{line.count(chr(9)) + 1} tab-separated fields, header has {len(fields)}. "
                    + _HINT
                )
            _check_no_nul(input_file, lineno, line)
            dst.write(f"{row}\t{line}\n")


def ensure_sorted_by_gene(input_file: Path, gene_col: str, context: PipelineContext) -> Path:
    """Ensure the TSV file is sorted by gene column for gene-aware chunking."""
    if context.config.get("assume_sorted", False):
        return input_file

    sorted_file = context.workspace.intermediate_dir / f"{input_file.stem}.sorted.tsv"

    if sorted_file.exists():
        logger.info(f"Using existing sorted file: {sorted_file}")
        return sorted_file

    logger.info(f"Sorting {input_file} by {gene_col} column")

    if str(input_file).endswith(".gz"):

        def open_func(f: Any) -> Any:
            return gzip.open(f, "rt")

    else:

        def open_func(f: Any) -> Any:
            return open(f)

    with open_func(input_file) as f:
        header = f.readline().strip().split("\t")
        try:
            col_idx = header.index(gene_col) + 1
        except ValueError:
            logger.error(f"Column {gene_col} not found in header")
            return input_file

    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, dir=context.config.get("temp_dir")
    ) as tmp:
        with open_func(input_file) as f:
            tmp.write(f.readline())
        tmp_path = tmp.name

    try:
        sort_cmd = [
            "sort",
            "-t",
            "\t",
            "-k",
            f"{col_idx},{col_idx}",  # the gene field only, not the rest of the line
            "--stable",
            "-T",
            str(context.workspace.temp_dir),
        ]

        file_size_mb = input_file.stat().st_size / (1024 * 1024)
        if file_size_mb > 1000:
            memory_limit = context.config.get("sort_memory_limit", "2G")
        elif file_size_mb > 100:
            memory_limit = context.config.get("sort_memory_limit", "1G")
        else:
            memory_limit = context.config.get("sort_memory_limit", "512M")

        if memory_limit:
            sort_cmd.extend(["-S", memory_limit])

        logger.info(f"Sorting {input_file} ({file_size_mb:.1f}MB) with {memory_limit} memory limit")

        with open(tmp_path, "a") as out_f:
            if str(input_file).endswith(".gz"):
                zcat_proc = subprocess.Popen(["zcat", str(input_file)], stdout=subprocess.PIPE)
                tail_proc = subprocess.Popen(
                    ["tail", "-n", "+2"], stdin=zcat_proc.stdout, stdout=subprocess.PIPE
                )
                assert zcat_proc.stdout is not None
                zcat_proc.stdout.close()
            else:
                tail_proc = subprocess.Popen(
                    ["tail", "-n", "+2", str(input_file)], stdout=subprocess.PIPE
                )

            assert tail_proc.stdout is not None
            # Bytewise collation: locale rules could interleave e.g. MYC and MYCN rows.
            sort_env = {**os.environ, "LC_ALL": "C"}
            sort_proc = subprocess.Popen(
                sort_cmd, stdin=tail_proc.stdout, stdout=out_f, env=sort_env
            )
            tail_proc.stdout.close()
            sort_proc.wait()

            if sort_proc.returncode != 0:
                raise subprocess.CalledProcessError(sort_proc.returncode, sort_cmd)

        shutil.move(tmp_path, sorted_file)
        logger.info(f"Sorting complete: {sorted_file}")
        return sorted_file

    except Exception as e:
        logger.error(f"Sorting failed: {e}")
        with contextlib.suppress(Exception):
            os.unlink(tmp_path)
        return input_file


def calculate_inheritance_safe_chunk_size(
    sample_count: int, context: PipelineContext | None = None
) -> int:
    """Calculate safe chunk size based on available memory to prevent inheritance skipping."""
    logger.debug(f"Calculating inheritance-safe chunk size for {sample_count} samples")

    try:
        if context and hasattr(context, "config"):
            rm = getattr(context, "resource_manager", None)
            if rm is None:
                from ..memory import ResourceManager

                rm = ResourceManager(config=context.config)
            total_variants = 100000
            chunk_size = rm.auto_chunk_size(total_variants, sample_count)
            logger.info(
                f"Using memory-aware chunk size: {chunk_size} variants for {sample_count} samples"
            )
            return chunk_size
    except Exception as e:
        logger.debug(f"Could not use ResourceManager: {e}. Falling back to improved defaults.")

    if sample_count <= 100:
        chunk_size = 50000
    elif sample_count <= 500:
        chunk_size = 25000
    elif sample_count <= 1000:
        chunk_size = 15000
    elif sample_count <= 2500:
        chunk_size = 8000
    elif sample_count <= 5000:
        chunk_size = 6000
    else:
        target_cells = 50_000_000
        chunk_size = max(1000, target_cells // sample_count)

    chunk_size = max(1000, min(chunk_size, 100000))

    logger.info(
        f"Using inheritance-safe chunk size: {chunk_size} variants for {sample_count} samples "
        f"(~{chunk_size * sample_count / 1_000_000:.1f}M cells)"
    )
    return chunk_size


def read_tsv_in_gene_chunks(
    filepath: Path,
    gene_column: str,
    chunksize: int,
    context: PipelineContext | None = None,
    max_variants_per_gene: int = 10000,
) -> Iterator[pd.DataFrame]:
    """Read TSV file in gene-aware chunks with configurable gene size limits."""
    if context and hasattr(context, "vcf_samples") and context.vcf_samples:
        sample_count = len(context.vcf_samples)
        safe_chunk_size = calculate_inheritance_safe_chunk_size(sample_count, context)

        if chunksize > safe_chunk_size:
            logger.warning(
                f"User-specified chunk size ({chunksize}) exceeds inheritance-safe limit "
                f"({safe_chunk_size}) for {sample_count} samples. Using safe chunk size."
            )
            chunksize = safe_chunk_size
        else:
            logger.info(
                f"Using user-specified chunk size ({chunksize}) which is within "
                f"inheritance-safe limit ({safe_chunk_size}) for {sample_count} samples"
            )
    else:
        logger.debug("No VCF samples available, using original chunk size")

    logger.info(
        f"Reading {filepath} in gene-aware chunks (column: {gene_column}, size: {chunksize})"
    )

    reader = pd.read_csv(filepath, chunksize=chunksize, **TSV_READ_KWARGS)

    gene_buffer = pd.DataFrame()
    chunks_yielded = 0

    if context and hasattr(context, "config"):
        max_variants_per_gene = context.config.get("max_variants_per_gene", max_variants_per_gene)

    max_buffer_size = min(chunksize * 5, max_variants_per_gene)
    emergency_threshold = min(chunksize * 10, max_variants_per_gene)

    for chunk in reader:
        if gene_buffer.empty:
            gene_buffer = chunk
        else:
            gene_buffer = pd.concat([gene_buffer, chunk], ignore_index=True)

        if len(gene_buffer) > max_buffer_size:
            logger.warning(
                f"Gene buffer has grown to {len(gene_buffer)} rows. "
                f"This may indicate a gene with too many variants."
            )
            if len(gene_buffer) > emergency_threshold:
                logger.error(
                    f"Gene buffer exceeded emergency threshold ({emergency_threshold} rows). "
                    f"Forcing chunk yield to prevent memory exhaustion."
                )
                gene_values = gene_buffer[gene_column].values
                for i in range(chunksize, len(gene_values)):
                    if gene_values[i] != gene_values[i - 1]:
                        yield_df = gene_buffer.iloc[:i].copy()
                        yield yield_df
                        chunks_yielded += 1
                        gene_buffer = gene_buffer.iloc[i:].reset_index(drop=True)
                        break
                else:
                    logger.warning("No gene boundaries found, yielding entire buffer")
                    yield gene_buffer.copy()
                    chunks_yielded += 1
                    gene_buffer = pd.DataFrame()

        if len(gene_buffer) >= chunksize:
            gene_values = gene_buffer[gene_column].values
            gene_change_indices = [0]
            current_gene = gene_values[0]

            for i in range(1, len(gene_values)):
                if gene_values[i] != current_gene:
                    gene_change_indices.append(i)
                    current_gene = gene_values[i]

            split_index = None
            for idx in gene_change_indices:
                if idx >= chunksize:
                    split_index = idx
                    break

            if split_index is None:
                continue
            else:
                yield_df = gene_buffer.iloc[:split_index].copy()
                if logger.isEnabledFor(logging.DEBUG):
                    genes_in_chunk = yield_df[gene_column].nunique()
                    logger.debug(
                        f"Yielding chunk with {len(yield_df)} variants across "
                        f"{genes_in_chunk} genes"
                    )

                yield yield_df
                chunks_yielded += 1
                gene_buffer = gene_buffer.iloc[split_index:].reset_index(drop=True)

    if not gene_buffer.empty:
        if logger.isEnabledFor(logging.DEBUG):
            genes_in_final_chunk = gene_buffer[gene_column].nunique()
            logger.debug(
                f"Yielding final chunk with {len(gene_buffer)} variants across "
                f"{genes_in_final_chunk} genes"
            )
        yield gene_buffer
        chunks_yielded += 1

    logger.info(f"Completed reading {chunks_yielded} gene-aware chunks")


def apply_chunk_annotations(chunk_df: pd.DataFrame, context: PipelineContext) -> pd.DataFrame:
    """Apply custom annotations to a chunk."""
    try:
        from ..annotator import annotate_dataframe_with_features, load_custom_features

        annotation_cfg = {
            "annotate_bed_files": context.annotation_configs.get("bed_files", []),
            "annotate_gene_lists": context.annotation_configs.get("gene_lists", []),
            "annotate_json_genes": context.annotation_configs.get("json_genes", []),
            "json_gene_mapping": context.annotation_configs.get("json_mapping", ""),
            "json_genes_as_columns": context.config.get("json_genes_as_columns", False),
            "vcf_samples": context.vcf_samples,
        }
        features = load_custom_features(annotation_cfg)

        if features:
            chunk_df = annotate_dataframe_with_features(chunk_df, features)
            logger.debug(f"Applied {len(features)} annotation features to chunk")

        return chunk_df
    except Exception as e:
        logger.error(f"Error applying annotations to chunk: {e}")
        raise


def apply_chunk_inheritance_analysis(
    chunk_df: pd.DataFrame, context: PipelineContext
) -> pd.DataFrame:
    """Apply inheritance analysis to a chunk with gene-aware processing."""
    try:
        from ..inheritance.analyzer import analyze_inheritance, process_inheritance_output

        inheritance_config = context.config.get("inheritance_analysis_config", {})
        vcf_samples = inheritance_config.get("vcf_samples", [])
        pedigree_data = inheritance_config.get("pedigree_data", {})
        use_vectorized_comp_het = inheritance_config.get("use_vectorized_comp_het", True)
        inheritance_mode = inheritance_config.get("inheritance_mode", "simple")
        preserve_details_for_scoring = inheritance_config.get("preserve_details_for_scoring", False)

        if not vcf_samples:
            logger.warning("No VCF samples available for chunk inheritance analysis")
            return chunk_df

        rm = context.resource_manager
        if rm is None:
            from ..memory import ResourceManager

            rm = ResourceManager(config=context.config)
        estimated_memory_gb = rm.estimate_memory(len(chunk_df), len(vcf_samples))
        safe_memory_gb = rm.memory_gb * rm.memory_safety_factor

        logger.debug(
            f"Chunk memory: {estimated_memory_gb:.2f}GB estimated, {safe_memory_gb:.2f}GB available"
        )

        force_processing = context.config.get("force_inheritance_processing", False)
        should_process = estimated_memory_gb <= safe_memory_gb or force_processing

        if not should_process:
            logger.error(
                f"MEMORY LIMIT: Skipping inheritance analysis for chunk with "
                f"{len(chunk_df)} variants x {len(vcf_samples)} samples "
                f"({estimated_memory_gb:.2f}GB exceeds {safe_memory_gb:.2f}GB limit)."
            )
            chunk_df["Inheritance_Pattern"] = "memory_limit_exceeded"
            chunk_df["Inheritance_Details"] = '{"error": "memory_limit_exceeded"}'
            chunk_df = process_inheritance_output(
                chunk_df,
                inheritance_mode,
                preserve_details_for_scoring=preserve_details_for_scoring,
            )
            return cleanup_sample_columns(
                chunk_df,
                vcf_samples,
                preserve_columns=["GT", "Inheritance_Pattern", "Inheritance_Details"],
            )

        logger.info(
            f"Starting inheritance analysis for chunk: {len(chunk_df)} variants, "
            f"{len(vcf_samples)} samples, inheritance mode: {inheritance_mode}"
        )

        sample_columns_exist = any(s in chunk_df.columns for s in vcf_samples)
        if not sample_columns_exist and "GT" in chunk_df.columns and len(chunk_df) > 0:
            logger.debug("Creating sample columns from GT column for inheritance analysis")
            prep_start = time.time()
            sample_column_method = context.config.get("sample_column_creation_method", "auto")
            chunk_df = create_sample_columns_from_gt_intelligent(
                df=chunk_df,
                vcf_samples=vcf_samples,
                separator=context.config.get("separator", ";"),
                snpsift_sep=context.config.get("extract_fields_separator", ","),
                method=sample_column_method,
            )
            prep_time = time.time() - prep_start
            if prep_time > 0.1:
                logger.debug(f"Sample column preparation took {prep_time:.2f}s")
        elif not sample_columns_exist:
            gt_cols = _find_per_sample_gt_columns(chunk_df)
            if gt_cols:
                n = _copy_gt_columns_to_samples(chunk_df, gt_cols, vcf_samples)  # by GEN index
                logger.info(f"Created sample columns from {n} per-sample GT columns")
            else:
                logger.warning("No GT column or per-sample GT columns found in chunk")

        analysis_start = time.time()
        threads = context.config.get("threads", 1)
        min_variants_for_parallel = context.config.get("min_variants_for_parallel_inheritance", 100)

        if threads > 1 and len(chunk_df) >= min_variants_for_parallel:
            logger.info(f"Using parallel inheritance analyzer for chunk with {threads} workers")
            from ..inheritance.parallel_analyzer import analyze_inheritance_parallel

            comp_het_start = time.time()
            chunk_df = analyze_inheritance_parallel(
                df=chunk_df,
                sample_list=vcf_samples,
                pedigree_data=pedigree_data,
                use_vectorized_comp_het=use_vectorized_comp_het,
                n_workers=threads,
                min_variants_for_parallel=min_variants_for_parallel,
                reference=context.config.get("reference"),
            )
            comp_het_time = time.time() - comp_het_start
            if comp_het_time > 1.0:
                logger.debug(f"Parallel compound het analysis took {comp_het_time:.2f}s")
        else:
            chunk_df = analyze_inheritance(
                df=chunk_df,
                pedigree_data=pedigree_data,
                sample_list=vcf_samples,
                use_vectorized_comp_het=use_vectorized_comp_het,
                reference=context.config.get("reference"),
            )

        analysis_time = time.time() - analysis_start
        logger.info(
            f"Inheritance analysis completed in {analysis_time:.2f}s for {len(chunk_df)} variants"
        )

        output_start = time.time()
        if preserve_details_for_scoring:
            logger.debug("Preserving Inheritance_Details for chunk scoring")

        chunk_df = process_inheritance_output(
            chunk_df,
            inheritance_mode,
            preserve_details_for_scoring=preserve_details_for_scoring,
        )

        output_time = time.time() - output_start
        if output_time > 0.1:
            logger.debug(f"Output processing took {output_time:.2f}s")

        cleanup_start = time.time()
        chunk_df = cleanup_sample_columns(
            chunk_df,
            vcf_samples,
            preserve_columns=["GT", "Inheritance_Pattern", "Inheritance_Details"],
        )

        cleanup_time = time.time() - cleanup_start
        if cleanup_time > 0.1:
            logger.debug(f"Column cleanup took {cleanup_time:.2f}s")

        inheritance_patterns = chunk_df["Inheritance_Pattern"].value_counts()
        logger.info(
            f"Chunk inheritance analysis SUCCESS: {len(chunk_df)} variants processed. "
            f"Patterns: {dict(inheritance_patterns)}"
        )
        return chunk_df

    except Exception as e:
        logger.error(f"Error in chunk inheritance analysis: {e}")
        raise

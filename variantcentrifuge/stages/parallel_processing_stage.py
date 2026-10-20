"""Parallel complete processing pipeline stage for variant chunks."""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
from pathlib import Path
from typing import Any

from ..extractor import extract_fields_bcftools, get_expected_extraction_columns
from ..filters import apply_snpsift_filter, extract_variants
from ..memory.resource_manager import ResourceManager
from ..pipeline_core import PipelineContext, Stage
from ..sample_subset import resolve_sample_subset
from ..transcript_filter import filter_vcf_to_transcripts, load_transcript_ids
from ..utils import ensure_fields_in_extract, split_bed_file
from ..vcf_eff_one_per_line import process_vcf_file as split_snpeff_annotations

logger = logging.getLogger(__name__)

__all__ = ["ParallelCompleteProcessingStage"]


def _ps_attr(name: str, fallback: Any) -> Any:
    """Get attribute from this module or processing_stages (supports test patching)."""
    ps_mod = sys.modules.get("variantcentrifuge.stages.processing_stages")
    if ps_mod is not None and hasattr(getattr(ps_mod, name, None), "_mock_name"):
        return getattr(ps_mod, name)
    this_mod = sys.modules.get(__name__)
    if this_mod is not None and getattr(this_mod, name, None) is not None:
        return getattr(this_mod, name)
    return fallback


class ParallelCompleteProcessingStage(Stage):
    """Run complete processing pipeline in parallel for each BED chunk."""

    @property
    def name(self) -> str:
        return "parallel_complete_processing"

    @property
    def checkpoint_resume_policy(self) -> str:
        return "restore"

    @property
    def description(self) -> str:
        return "Process variants in parallel with complete pipeline per chunk"

    @property
    def dependencies(self) -> set[str]:
        return {"gene_bed_creation", "configuration_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        # Sample/phenotype stages fix the analysed sample set used by resolve_sample_subset.
        return {"sample_config_loading", "phenotype_case_control_assignment"}

    @property
    def parallel_safe(self) -> bool:
        return False

    @property
    def estimated_runtime(self) -> float:
        return 60.0

    def _chunk_fingerprint(self, context: PipelineContext, chunk_bed: Path, eff_fields: str) -> str:
        """Compute a short fingerprint of chunk BED content and extraction/filtering config."""
        cfg = context.config
        bed_bytes = chunk_bed.read_bytes() if chunk_bed.is_file() else str(chunk_bed).encode()
        keys = (
            "bcftools_prefilter",
            "late_filtering",
            "snpeff_splitting_mode",
            "split_snpeff_lines",
            "filter",
            "filters",
            "filter_exception_file",
            "gene_name",
            "gene_file",
            "regions_bed",
            "transcript_list",
            "transcript_file",
        )
        data = {k: cfg.get(k) for k in keys}
        data["bed"] = hashlib.sha256(bed_bytes).hexdigest()[:16]
        data["vcf"] = str(cfg.get("vcf_file") or context.data or "")
        data["fields"] = eff_fields
        data["samples"] = cfg.get("sample_subset_sha256")
        return hashlib.sha256(json.dumps(data, sort_keys=True, default=str).encode()).hexdigest()[
            :16
        ]

    def _validate_chunk_tsv(
        self,
        chunk_tsv: Path,
        fallback_on_error: bool = True,
        expected_header: str | None = None,
        expected_fingerprint: str | None = None,
    ) -> bool:
        """Validate that a chunk TSV file is complete and valid."""
        if not chunk_tsv.exists():
            return False
        if expected_fingerprint is not None:
            meta = Path(f"{chunk_tsv}.meta")
            if (
                not meta.is_file()
                or meta.read_text(encoding="utf-8").strip() != expected_fingerprint
            ):
                return False
        try:
            if chunk_tsv.stat().st_size == 0:
                return False
            open_func = gzip.open if chunk_tsv.name.endswith(".gz") else open
            mode = "rt" if chunk_tsv.name.endswith(".gz") else "r"
            with open_func(chunk_tsv, mode) as f:
                first_line = str(f.readline()).strip()
                if not first_line or "\t" not in first_line:
                    return False
                if expected_header is not None and first_line != expected_header.strip():
                    return False
            return True
        except Exception:
            if not fallback_on_error:
                raise
            return False

    def _build_expected_header(self, context: PipelineContext) -> str | None:
        """Build expected TSV header string from context config and samples."""
        eff_fields = self._resolve_effective_extraction_fields(context.config)
        if not eff_fields:
            return None
        get_cols = _ps_attr("get_expected_extraction_columns", get_expected_extraction_columns)
        return "\t".join(get_cols(eff_fields, getattr(context, "vcf_samples", None)))

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system."""
        resolve_sample_subset(context)  # the restored TSV holds only the subset's columns
        ws = context.workspace
        gz_tsv = ws.intermediate_dir / f"{ws.base_name}.extracted.tsv.gz"
        tsv = ws.intermediate_dir / f"{ws.base_name}.extracted.tsv"
        merged_tsv = gz_tsv if gz_tsv.exists() else tsv
        exp_header = self._build_expected_header(context)
        if not merged_tsv.exists() or not self._validate_chunk_tsv(
            merged_tsv, expected_header=exp_header
        ):
            raise RuntimeError(
                f"Cannot restore parallel_complete_processing from checkpoint: "
                f"merged TSV missing or invalid: {merged_tsv}"
            )
        context.extracted_tsv = merged_tsv
        context.data = merged_tsv
        context.config["parallel_vcf_processing_complete"] = True
        context.config.pop("chunked_processing_complete", None)
        for stage in ("variant_extraction", "snpsift_filtering", "field_extraction"):
            context.mark_complete(stage)
        logger.info(f"Restored valid merged TSV: {merged_tsv}")
        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process variants in parallel with complete pipeline per chunk."""
        validate_fn = _ps_attr("validate_annotation_consequence_filtering", None)
        if validate_fn is None:
            from .processing_stages import validate_annotation_consequence_filtering as validate_fn
        validate_fn(context.config)
        threads = context.config.get("threads", 1)

        if threads <= 1:
            logger.info("Thread count <= 1, falling back to sequential stages")
            return context

        if context.gene_bed_file is None:
            raise ValueError("gene_bed_file is required for parallel complete processing")
        resolve_sample_subset(context)  # before chunk headers/fingerprints use vcf_samples
        split_start = self._start_subtask("bed_splitting")
        # Chunks live in the run's workspace: a user --regions-bed may sit in a read-only or
        # shared input directory, and concurrent runs on one BED must not collide.
        bed_chunks = self._split_bed_file(
            context.gene_bed_file, threads, context.workspace.intermediate_dir / "chunks"
        )
        self._end_subtask("bed_splitting", split_start)
        logger.info(f"Split BED file into {len(bed_chunks)} chunks for parallel processing")

        process_start = self._start_subtask("parallel_chunk_processing")
        chunk_tsvs = self._process_chunks_parallel(context, bed_chunks)
        self._end_subtask("parallel_chunk_processing", process_start)

        merge_start = self._start_subtask("tsv_merging")
        merged_tsv = self._merge_tsv_outputs(context, chunk_tsvs)
        self._end_subtask("tsv_merging", merge_start)

        context.extracted_tsv = merged_tsv
        context.data = merged_tsv
        for stage in ("variant_extraction", "snpsift_filtering", "field_extraction"):
            context.mark_complete(stage)
        context.config["parallel_vcf_processing_complete"] = True
        context.config.pop("chunked_processing_complete", None)

        cleanup_start = self._start_subtask("cleanup")
        self._cleanup_chunks(bed_chunks, chunk_tsvs, context)
        self._end_subtask("cleanup", cleanup_start)
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        files = [context.gene_bed_file] if context.gene_bed_file is not None else []
        vcf_file = getattr(context, "vcf_file", None)
        if vcf_file:
            files.append(Path(vcf_file))
        return files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        return [context.extracted_tsv] if context.extracted_tsv is not None else []

    def _split_bed_file(self, bed_file: Path, n_chunks: int, chunk_dir: Path) -> list[Path]:
        """Split BED file into roughly equal chunks inside ``chunk_dir``."""
        chunk_dir.mkdir(exist_ok=True)
        split_fn = _ps_attr("split_bed_file", split_bed_file)
        chunks = split_fn(str(bed_file), n_chunks, str(chunk_dir))
        return [Path(chunk) for chunk in chunks]

    def _resolve_effective_extraction_fields(self, config: dict[str, Any]) -> str:
        """Resolve effective extraction field string with extra sample fields and RAW_FILTER."""
        fields_str = config.get("fields_to_extract", "")
        ensure_fn = _ps_attr("ensure_fields_in_extract", ensure_fields_in_extract)
        if config.get("append_extra_sample_fields", False) and config.get("extra_sample_fields"):
            fields_str = ensure_fn(fields_str, config["extra_sample_fields"])
        if config.get("filter_exception_file"):
            fields_str = ensure_fn(fields_str, ["RAW_FILTER"])
        return str(fields_str)

    def _calculate_chunk_worker_resources(
        self,
        total_threads: int,
        num_chunks: int,
        shared_memory_gb: float,
        requires_snpsift: bool = True,
    ) -> dict[str, Any]:
        """Calculate active worker count, threads per worker, and SnpSift heap per worker."""
        if num_chunks <= 0:
            return {
                "active_workers": 0,
                "threads_per_worker": max(1, total_threads),
                "snpsift_memory_gb": 2,
                "max_memory_gb": shared_memory_gb,
            }
        if not requires_snpsift:
            active_workers = max(1, min(num_chunks, max(1, total_threads)))
            return {
                "active_workers": active_workers,
                "threads_per_worker": max(1, total_threads // active_workers),
                "snpsift_memory_gb": 2,
                "max_memory_gb": max(0.5, shared_memory_gb / active_workers),
            }
        total_snpsift_gb = int(shared_memory_gb * 0.5)
        if total_snpsift_gb < 2:
            raise ValueError(
                f"Memory budget ({shared_memory_gb} GiB) is too small to accommodate the "
                "2 GiB minimum SnpSift heap plus non-heap headroom (at least 4 GiB required)."
            )
        desired_workers = min(num_chunks, max(1, total_threads))
        active_workers = max(1, min(desired_workers, max(1, total_snpsift_gb // 2)))
        return {
            "active_workers": active_workers,
            "threads_per_worker": max(1, total_threads // active_workers),
            "snpsift_memory_gb": max(2, min(16, total_snpsift_gb // active_workers)),
            "max_memory_gb": max(1.0, shared_memory_gb / active_workers),
        }

    def _process_single_chunk(
        self,
        chunk_index: int,
        chunk_bed: Path,
        vcf_file: str,
        base_name: str,
        intermediate_dir: Path,
        config: dict[str, Any],
        subtask_times: Any = None,
    ) -> Path:
        """Process a single BED chunk through complete pipeline."""
        chunk_base = f"{base_name}.chunk_{chunk_index}"
        extract_start = time.time()
        chunk_vcf = intermediate_dir / f"{chunk_base}.variants.vcf.gz"
        _ps_attr("extract_variants", extract_variants)(
            vcf_file=vcf_file,
            bed_file=str(chunk_bed),
            cfg={
                "threads": config.get("threads_per_chunk", 1),
                "bcftools_prefilter": config.get("bcftools_prefilter"),
            },
            output_file=str(chunk_vcf),
        )
        extract_time = time.time() - extract_start

        working_vcf = chunk_vcf
        split_mode = config.get("snpeff_splitting_mode")
        split_fn = _ps_attr("split_snpeff_annotations", split_snpeff_annotations)
        if split_mode and split_mode != "after_filters":
            split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
            split_fn(str(working_vcf), str(split_vcf))
            working_vcf = split_vcf

        transcript_ids = config.get("transcript_ids")
        if transcript_ids is None:
            load_tx_fn = _ps_attr("load_transcript_ids", load_transcript_ids)
            transcript_ids = load_tx_fn(
                config.get("transcript_list"), config.get("transcript_file")
            )
        if transcript_ids:
            transcript_filtered_vcf = intermediate_dir / f"{chunk_base}.transcript_filtered.vcf.gz"
            _ps_attr("filter_vcf_to_transcripts", filter_vcf_to_transcripts)(
                str(working_vcf),
                str(transcript_filtered_vcf),
                transcript_ids,
                index_output=False,
                filter_exception_file=config.get("filter_exception_file"),
            )
            working_vcf = transcript_filtered_vcf

        filter_start = time.time()
        if config.get("late_filtering", False):
            filtered_vcf = working_vcf
            filter_time = 0.0
        else:
            filtered_vcf = intermediate_dir / f"{chunk_base}.filtered.vcf.gz"
            filter_expr = config.get("filter") or config.get("filters")
            if filter_expr:
                _ps_attr("apply_snpsift_filter", apply_snpsift_filter)(
                    str(working_vcf),
                    filter_expr,
                    {
                        "threads": config.get("threads_per_chunk", 1),
                        "max_memory_gb": config.get("max_memory_gb"),
                        "snpsift_memory_gb": config.get("snpsift_memory_gb"),
                        "filter_exception_file": config.get("filter_exception_file"),
                    },
                    str(filtered_vcf),
                )
            else:
                filtered_vcf = working_vcf
        if split_mode == "after_filters":
            split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
            split_fn(str(filtered_vcf), str(split_vcf))
            filtered_vcf = split_vcf
        filter_time = time.time() - filter_start

        field_extract_start = time.time()
        use_compression = config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"
        chunk_tsv = intermediate_dir / f"{chunk_base}{tsv_suffix}"

        fields_str = self._resolve_effective_extraction_fields(config)
        fields = fields_str.split() if fields_str else []
        if not fields:
            raise ValueError("No fields specified for extraction")

        _ps_attr("extract_fields_bcftools", extract_fields_bcftools)(
            variant_file=str(filtered_vcf),
            fields=" ".join(fields),
            cfg={
                "threads": config.get("threads_per_chunk", 1),
                "extract_fields_separator": config.get("extract_fields_separator", ","),
                "bcftools_samples": config.get("bcftools_samples"),
            },
            output_file=str(chunk_tsv),
            vcf_samples=config.get("vcf_samples"),
        )
        field_extract_time = time.time() - field_extract_start

        if subtask_times is not None:
            subtask_times[f"chunk_{chunk_index}_extraction"] = extract_time
            subtask_times[f"chunk_{chunk_index}_filtering"] = filter_time
            subtask_times[f"chunk_{chunk_index}_field_extraction"] = field_extract_time
        logger.debug(f"Completed processing chunk {chunk_index}")
        return chunk_tsv

    def _process_chunks_parallel(
        self, context: PipelineContext, bed_chunks: list[Path]
    ) -> list[Path]:
        """Process BED chunks in parallel with complete pipeline and automatic substep detection."""
        if not bed_chunks:
            logger.warning("No BED chunks to process - skipping parallel processing")
            return []
        vcf_file = context.config.get("vcf_file") or (str(context.data) if context.data else None)
        if not vcf_file:
            raise KeyError("vcf_file")
        threads = context.config.get("threads", 1)
        base_name = context.workspace.base_name
        intermediate_dir = context.workspace.intermediate_dir
        chunks_to_process = []
        use_compression = context.config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"

        expected_header: str | None = None
        eff_fields = self._resolve_effective_extraction_fields(context.config)
        if eff_fields:
            get_cols_fn = _ps_attr(
                "get_expected_extraction_columns", get_expected_extraction_columns
            )
            cols = get_cols_fn(eff_fields, getattr(context, "vcf_samples", None))
            expected_header = "\t".join(cols)

        for i, chunk_bed in enumerate(bed_chunks):
            expected_tsv = intermediate_dir / f"{base_name}.chunk_{i}{tsv_suffix}"
            chunk_fp = self._chunk_fingerprint(context, chunk_bed, eff_fields)
            if self._validate_chunk_tsv(
                expected_tsv, expected_header=expected_header, expected_fingerprint=chunk_fp
            ):
                logger.info(f"Reusing existing chunk TSV {i}: {expected_tsv.name}")
            else:
                logger.info(f"Processing missing/invalid chunk TSV {i}")
                chunks_to_process.append((i, chunk_bed, expected_tsv, chunk_fp))

        if chunks_to_process:
            cfg = context.config
            rm = ResourceManager(config=cfg)
            shared_memory_gb = float(cfg.get("max_memory_gb") or rm.memory_gb)
            req_snpsift = bool(cfg.get("filter") or cfg.get("filters")) and not bool(
                cfg.get("late_filtering")
            )
            resource_alloc = self._calculate_chunk_worker_resources(
                total_threads=threads,
                num_chunks=len(chunks_to_process),
                shared_memory_gb=shared_memory_gb,
                requires_snpsift=req_snpsift,
            )
            active_workers = resource_alloc["active_workers"]
            threads_per_worker = resource_alloc["threads_per_worker"]
            worker_snpsift_gb = resource_alloc["snpsift_memory_gb"]
            worker_max_memory_gb = resource_alloc["max_memory_gb"]
            logger.info(
                f"Processing {len(chunks_to_process)} missing chunks with "
                f"{active_workers} parallel workers "
                f"({threads_per_worker} threads/worker, -Xmx{worker_snpsift_gb}g heap/worker)"
            )
            load_tx_fn = _ps_attr("load_transcript_ids", load_transcript_ids)
            worker_config = {
                "threads_per_chunk": threads_per_worker,
                "max_memory_gb": worker_max_memory_gb,
                "snpsift_memory_gb": worker_snpsift_gb,
                "bcftools_prefilter": cfg.get("bcftools_prefilter"),
                "late_filtering": cfg.get("late_filtering", False),
                "snpeff_splitting_mode": cfg.get("snpeff_splitting_mode"),
                "filter": cfg.get("filter"),
                "filters": cfg.get("filters"),
                "filter_exception_file": cfg.get("filter_exception_file"),
                "fields_to_extract": eff_fields,
                "append_extra_sample_fields": bool(cfg.get("append_extra_sample_fields")),
                "extra_sample_fields": cfg.get("extra_sample_fields", []),
                "extract_fields_separator": cfg.get("extract_fields_separator", ","),
                "gzip_intermediates": use_compression,
                "vcf_samples": getattr(context, "vcf_samples", None),
                "bcftools_samples": cfg.get("_bcftools_sample_subset"),
                "transcript_ids": load_tx_fn(
                    cfg.get("transcript_list"), cfg.get("transcript_file")
                ),
            }
            shared_subtask_times = Manager().dict()
            executor_cls = _ps_attr("ProcessPoolExecutor", ProcessPoolExecutor)
            as_completed_fn = _ps_attr("as_completed", as_completed)
            with executor_cls(max_workers=active_workers) as executor:
                future_to_chunk = {
                    executor.submit(
                        self._process_single_chunk,
                        i,
                        chunk_bed,
                        vcf_file,
                        base_name,
                        intermediate_dir,
                        worker_config,
                        shared_subtask_times,
                    ): (i, chunk_bed, expected_tsv, chunk_fp)
                    for i, chunk_bed, expected_tsv, chunk_fp in chunks_to_process
                }
                for future in as_completed_fn(future_to_chunk):
                    chunk_index, chunk_bed, _expected_tsv, chunk_fp = future_to_chunk[future]
                    try:
                        chunk_tsv = future.result()
                        if chunk_tsv.parent.exists():
                            Path(f"{chunk_tsv}.meta").write_text(chunk_fp, encoding="utf-8")
                        logger.debug(f"Completed chunk {chunk_index}: {chunk_bed.name}")
                    except Exception as e:
                        logger.error(f"Failed to process chunk {chunk_bed}: {e}")
                        raise
            if shared_subtask_times:
                self._aggregate_subtask_times(shared_subtask_times)
        else:
            logger.info("All chunk TSVs already exist and are valid - no processing needed")

        return [
            intermediate_dir / f"{base_name}.chunk_{i}{tsv_suffix}" for i in range(len(bed_chunks))
        ]

    def _aggregate_subtask_times(self, shared_subtask_times: Any) -> None:
        """Aggregate subtask times from parallel processing."""
        groups: dict[str, list[float]] = {
            "variant_extraction": [],
            "snpsift_filtering": [],
            "field_extraction": [],
        }
        for k, v in shared_subtask_times.items():
            if "_field_extraction" in k:
                groups["field_extraction"].append(v)
            elif "_extraction" in k:
                groups["variant_extraction"].append(v)
            elif "_filtering" in k:
                groups["snpsift_filtering"].append(v)
        for stage_name, times in groups.items():
            if times:
                self._subtask_times[stage_name] = sum(times) / len(times)

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

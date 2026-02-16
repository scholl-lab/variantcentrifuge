"""
Processing stages for variant extraction and data transformation.

This module contains stages that handle the core data processing tasks:
- Gene BED file creation
- Variant extraction from VCF files
- Filtering with bcftools and SnpSift
- Field extraction to TSV format
- Phenotype integration
"""

import contextlib
import gzip
import logging
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from smart_open import smart_open

from ..extractor import extract_fields_bcftools
from ..filters import apply_snpsift_filter, extract_variants
from ..gene_bed import get_gene_bed, normalize_genes
from ..phenotype import extract_phenotypes_for_gt_row, extract_phenotypes_from_sample_columns
from ..pipeline_core import PipelineContext, Stage
from ..pipeline_core.error_handling import (
    FileFormatError,
    ToolNotFoundError,
    graceful_error_handling,
    retry_on_failure,
    validate_file_exists,
)
from ..utils import ensure_fields_in_extract, run_command, split_bed_file
from ..vcf_eff_one_per_line import process_vcf_file as split_snpeff_annotations

logger = logging.getLogger(__name__)


class GeneBedCreationStage(Stage):
    """Convert gene names to BED file with genomic intervals."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "gene_bed_creation"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Create BED file from gene names"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - creates independent BED file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Create BED file for specified genes."""
        with graceful_error_handling(self.name):
            # Get normalized gene names
            gene_name = context.config.get("gene_name", "")
            gene_file = context.config.get("gene_file")

            # Validate gene file if provided
            if gene_file:
                gene_file = str(validate_file_exists(gene_file, self.name))

            # Normalize genes
            try:
                normalized_genes = normalize_genes(gene_name, gene_file, logger)
                context.config["normalized_genes"] = normalized_genes
            except SystemExit as e:
                # Convert SystemExit to proper exception
                # Check the specific error condition
                if not gene_name and not gene_file:
                    raise ValueError(
                        "No gene name provided. Use --gene-name or --gene-file to specify genes."
                    ) from e
                elif gene_file and not os.path.exists(gene_file):
                    raise FileNotFoundError(f"Gene file not found: {gene_file}") from e
                elif gene_name and os.path.exists(gene_name):
                    raise ValueError(
                        f"File path '{gene_name}' provided to --gene-name. "
                        f"Use --gene-file for gene list files."
                    ) from e
                else:
                    # Generic error for other SystemExit cases
                    raise ValueError("Invalid gene specification") from e

            # Generate BED file
            reference = context.config.get("reference", "GRCh37.75")
            interval_expand = context.config.get("interval_expand", 0)
            add_chr = context.config.get("add_chr", True)

            logger.info(f"Creating BED file for genes: {normalized_genes}")

            # Retry logic for external tool failures
            @retry_on_failure(
                max_attempts=3, delay=2.0, exceptions=(subprocess.CalledProcessError,)
            )
            def create_bed():
                return get_gene_bed(
                    reference=reference,
                    gene_name=normalized_genes,
                    interval_expand=interval_expand,
                    add_chr=add_chr,
                    output_dir=str(context.workspace.output_dir),
                )

            try:
                bed_file = create_bed()
            except subprocess.CalledProcessError as e:
                # Check if it's a tool not found error
                if "snpEff" in str(e) and ("not found" in str(e) or "No such file" in str(e)):
                    raise ToolNotFoundError("snpEff", self.name) from e
                # Check if it's a gene not found error
                if "not found in database" in str(e):
                    raise FileFormatError(
                        normalized_genes, f"valid gene names for {reference}", self.name
                    ) from e
                raise

            context.gene_bed_file = Path(bed_file)
            logger.info(f"Created BED file: {bed_file}")

            return context

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the generated BED file."""
        if context.gene_bed_file:
            return [context.gene_bed_file]
        return []


class VariantExtractionStage(Stage):
    """Extract variants from VCF file using BED regions (single-threaded)."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "variant_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract variants from VCF by genomic regions"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"gene_bed_creation"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract variants using the BED file."""
        with graceful_error_handling(self.name):
            # Validate inputs
            vcf_file = context.config.get("vcf_file")
            if not vcf_file:
                raise ValueError("No VCF file specified")

            vcf_path = validate_file_exists(vcf_file, self.name)

            if not context.gene_bed_file:
                raise ValueError("No BED file available from gene extraction")

            bed_file = str(context.gene_bed_file)

            # Output path
            output_vcf = context.workspace.get_intermediate_path(
                f"{context.workspace.base_name}.variants.vcf.gz"
            )

            # Prepare config for extract_variants
            extract_config = {
                "threads": context.config.get("threads", 1),
                "bcftools_prefilter": context.config.get("bcftools_prefilter"),
            }

            logger.info(f"Extracting variants from {vcf_file} using {bed_file}")

            @retry_on_failure(
                max_attempts=3, delay=2.0, exceptions=(subprocess.CalledProcessError,)
            )
            def extract():
                extract_variants(
                    vcf_file=str(vcf_path),
                    bed_file=bed_file,
                    cfg=extract_config,
                    output_file=str(output_vcf),
                )

            try:
                extract()
            except subprocess.CalledProcessError as e:
                # Check for specific error conditions
                if "bcftools" in str(e) and ("not found" in str(e) or "No such file" in str(e)):
                    raise ToolNotFoundError("bcftools", self.name) from e
                if "invalid" in str(e).lower() and "vcf" in str(e).lower():
                    raise FileFormatError(str(vcf_path), "valid VCF format", self.name) from e
                raise

            # Verify output was created
            if not output_vcf.exists():
                raise FileNotFoundError(f"Variant extraction failed to create output: {output_vcf}")

            context.extracted_vcf = output_vcf
            context.data = output_vcf

            return context

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the extracted VCF file."""
        if context.extracted_vcf:
            return [context.extracted_vcf]
        return []


class ParallelVariantExtractionStage(Stage):
    """Extract variants in parallel by splitting BED file into chunks."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_variant_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract variants in parallel using multiple genomic regions"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"gene_bed_creation"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        # This stage manages its own parallelism internally
        return False

    @property
    def estimated_runtime(self) -> float:
        """Return the estimated runtime in seconds."""
        return 30.0  # Typically takes longer

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the expected output
        and mark dependent stages as complete.
        """
        # Reconstruct the expected merged VCF file path
        expected_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.variants.vcf.gz"
        )

        if expected_vcf.exists():
            logger.info(f"Restored parallel variant extraction output: {expected_vcf}")
            context.extracted_vcf = expected_vcf
            context.data = expected_vcf

            # Mark variant_extraction as complete for dependent stages
            context.mark_complete("variant_extraction")
        else:
            logger.warning(f"Expected VCF file not found during checkpoint skip: {expected_vcf}")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract variants in parallel by processing BED chunks."""
        threads = context.config.get("threads", 1)

        if threads <= 1:
            # Fall back to single-threaded extraction
            logger.info("Thread count <= 1, using single-threaded extraction")
            extractor = VariantExtractionStage()
            return extractor._process(context)

        # Split BED file into chunks
        if context.gene_bed_file is None:
            raise ValueError("gene_bed_file is required for parallel variant extraction")
        bed_chunks = self._split_bed_file(context.gene_bed_file, threads)
        logger.info(f"Split BED file into {len(bed_chunks)} chunks for parallel processing")

        # Process chunks in parallel (with automatic substep detection)
        chunk_outputs = self._process_chunks_parallel(context, bed_chunks)

        # Merge results
        merged_output = self._merge_chunk_outputs(context, chunk_outputs)

        context.extracted_vcf = merged_output
        context.data = merged_output

        # Cleanup chunks
        self._cleanup_chunks(bed_chunks, chunk_outputs)

        # Mark variant_extraction as complete so field_extraction can proceed
        context.mark_complete("variant_extraction")

        return context

    def _split_bed_file(self, bed_file: Path, n_chunks: int) -> list[Path]:
        """Split BED file into roughly equal chunks."""
        chunk_dir = bed_file.parent / "chunks"
        chunk_dir.mkdir(exist_ok=True)

        chunks = split_bed_file(str(bed_file), n_chunks, str(chunk_dir))
        return [Path(chunk) for chunk in chunks]

    def _validate_existing_chunk(self, chunk_path: Path, fallback_on_error: bool = True) -> bool:
        """Validate that an existing chunk file is complete and valid.

        Parameters
        ----------
        chunk_path : Path
            Path to the chunk file to validate
        fallback_on_error : bool
            If True, validation errors result in False (file will be reprocessed)
            If False, validation errors raise exceptions

        Returns
        -------
        bool
            True if file is valid, False if invalid or missing
        """
        if not chunk_path.exists():
            logger.debug(f"Chunk {chunk_path} does not exist")
            return False

        try:
            # Size check (must be > 0)
            stat = chunk_path.stat()
            if stat.st_size == 0:
                logger.debug(f"Chunk {chunk_path} is empty, invalid")
                return False

            # For VCF files, try basic validation
            if chunk_path.suffix == ".gz" and chunk_path.name.endswith(".vcf.gz"):
                try:
                    # Quick check - try to read the first few lines
                    import gzip

                    with gzip.open(chunk_path, "rt") as f:
                        first_line = f.readline()
                        if not first_line.startswith("##fileformat=VCF"):
                            logger.debug(f"Chunk {chunk_path} has invalid VCF header")
                            return False
                except Exception as e:
                    logger.debug(f"Chunk {chunk_path} failed VCF header validation: {e}")
                    if not fallback_on_error:
                        raise
                    return False

            logger.debug(f"Chunk {chunk_path} is valid (size: {stat.st_size} bytes)")
            return True

        except Exception as e:
            logger.warning(f"Error validating chunk {chunk_path}: {e}")
            if not fallback_on_error:
                raise
            return False

    def _process_chunks_parallel(
        self, context: PipelineContext, bed_chunks: list[Path]
    ) -> list[Path]:
        """Process BED chunks in parallel with automatic substep detection."""
        vcf_file = context.config["vcf_file"]
        threads = context.config.get("threads", 1)

        # Phase 1: Check for existing valid chunks
        chunks_to_process = []
        existing_outputs = []

        for i, chunk_bed in enumerate(bed_chunks):
            expected_output = context.workspace.get_temp_path(f"chunk_{i}.variants.vcf.gz")

            if self._validate_existing_chunk(expected_output):
                logger.info(f"Reusing existing chunk {i}: {expected_output.name}")
                existing_outputs.append(expected_output)
            else:
                logger.info(f"Processing missing/invalid chunk {i}")
                chunks_to_process.append((i, chunk_bed, expected_output))

        # Phase 2: Process only missing chunks
        new_outputs = []
        if chunks_to_process:
            # Each worker gets limited threads
            threads_per_worker = max(1, threads // len(chunks_to_process))

            # Prepare config for extract_variants
            extract_config = {
                "threads": threads_per_worker,
                "bcftools_prefilter": context.config.get("bcftools_prefilter"),
            }

            logger.info(
                f"Processing {len(chunks_to_process)} missing chunks out of {len(bed_chunks)} total"
            )

            with ProcessPoolExecutor(max_workers=len(chunks_to_process)) as executor:
                # Submit jobs for missing chunks only
                future_to_chunk = {}
                for i, chunk_bed, output_vcf in chunks_to_process:
                    future = executor.submit(
                        extract_variants,
                        vcf_file=vcf_file,
                        bed_file=str(chunk_bed),
                        cfg=extract_config,
                        output_file=str(output_vcf),
                    )
                    future_to_chunk[future] = (i, chunk_bed, output_vcf)

                # Collect results
                for future in as_completed(future_to_chunk):
                    i, chunk_bed, output_vcf = future_to_chunk[future]
                    try:
                        future.result()  # Raises any exceptions
                        new_outputs.append(output_vcf)
                        logger.debug(f"Completed extraction for chunk {i}: {chunk_bed.name}")
                    except Exception as e:
                        logger.error(f"Failed to process chunk {chunk_bed}: {e}")
                        raise
        else:
            logger.info("All chunks already exist and are valid - no processing needed")

        # Phase 3: Combine existing and new outputs, maintaining order
        all_outputs = []
        for i in range(len(bed_chunks)):
            expected_output = context.workspace.get_temp_path(f"chunk_{i}.variants.vcf.gz")
            all_outputs.append(expected_output)

        return all_outputs

    def _merge_chunk_outputs(self, context: PipelineContext, chunk_outputs: list[Path]) -> Path:
        """Merge VCF chunks into a single file."""
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.variants.vcf.gz"
        )

        if len(chunk_outputs) == 1:
            # Just move the single file
            shutil.move(str(chunk_outputs[0]), str(output_vcf))
        else:
            # Merge with bcftools
            cmd = [
                "bcftools",
                "concat",
                "-a",  # Allow overlaps
                "-Oz",  # Output compressed
                "-o",
                str(output_vcf),
            ] + [str(f) for f in chunk_outputs]

            run_command(cmd)

            # Index the merged file
            run_command(["bcftools", "index", str(output_vcf)])

        logger.info(f"Merged {len(chunk_outputs)} chunks into {output_vcf}")
        return output_vcf

    def _cleanup_chunks(self, bed_chunks: list[Path], vcf_chunks: list[Path]) -> None:
        """Clean up temporary chunk files."""
        for chunk in bed_chunks + vcf_chunks:
            if chunk.exists():
                chunk.unlink()

        # Remove chunks directory if empty
        if bed_chunks and bed_chunks[0].parent.name == "chunks":
            chunk_dir = bed_chunks[0].parent
            if chunk_dir.exists() and not any(chunk_dir.iterdir()):
                chunk_dir.rmdir()


class BCFToolsPrefilterStage(Stage):
    """Apply bcftools pre-filtering for performance optimization.

    Note: This stage is typically not needed as pre-filtering is applied
    during the variant extraction stage for better performance.
    This remains here for cases where separate filtering is needed.
    """

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "bcftools_prefilter"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply bcftools filter expression"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"variant_extraction", "parallel_variant_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply bcftools filter if specified."""
        bcftools_prefilter = context.config.get("bcftools_prefilter")

        if not bcftools_prefilter:
            logger.debug("No bcftools prefilter specified, skipping")
            return context

        # Note: This is typically applied during extraction for efficiency
        # This stage is here for cases where we need a separate filtering step
        logger.info(f"Applying bcftools prefilter: {bcftools_prefilter}")

        input_vcf = context.extracted_vcf or context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.bcftools_filtered.vcf.gz"
        )

        cmd = [
            "bcftools",
            "view",
            "-i",
            bcftools_prefilter,
            "-Oz",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ]

        run_command(cmd)
        run_command(["bcftools", "index", str(output_vcf)])

        context.bcftools_filtered_vcf = output_vcf  # type: ignore[attr-defined]
        context.filtered_vcf = output_vcf
        # Only update context.data if we haven't extracted to TSV yet
        if not hasattr(context, "extracted_tsv"):
            context.data = output_vcf

        return context


class MultiAllelicSplitStage(Stage):
    """Split multi-allelic SNPeff annotations into one per line."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "multiallelic_split"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Split multi-allelic SNPeff annotations"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after variant extraction but before field extraction
        return {"variant_extraction"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - independent transformation, creates new file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Split SNPeff annotations if requested."""
        if not context.config.get("snpeff_split_by_transcript"):
            logger.debug("SNPeff transcript splitting not requested")
            return context

        # Determine when to split
        split_before = context.config.get("snpeff_split_before_filter", True)
        if split_before and context.is_complete("snpsift_filtering"):
            # Already filtered, don't split again
            return context

        input_vcf = context.filtered_vcf or context.extracted_vcf or context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.split_annotations.vcf.gz"
        )

        logger.info("Splitting SNPeff annotations by transcript")
        split_snpeff_annotations(str(input_vcf), str(output_vcf))

        context.split_annotations_vcf = output_vcf  # type: ignore[attr-defined]
        # Only update context.data if we haven't extracted to TSV yet
        if not hasattr(context, "extracted_tsv"):
            context.data = output_vcf
        return context


class TranscriptFilterStage(Stage):
    """Filter variants by transcript IDs using SnpSift."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "transcript_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Filter variants by transcript IDs"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after multiallelic split to ensure proper annotation structure
        return {"multiallelic_split", "variant_extraction"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Run after bcftools prefilter if present
        return {"bcftools_prefilter"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True  # Safe - independent transformation, creates new file

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply transcript filtering if requested."""
        # Parse transcripts from config
        transcripts = []

        # Get transcript list from config
        transcript_list = context.config.get("transcript_list")
        if transcript_list:
            transcripts.extend([t.strip() for t in transcript_list.split(",") if t.strip()])

        # Get transcript file from config
        transcript_file = context.config.get("transcript_file")
        if transcript_file:
            if not os.path.exists(transcript_file):
                logger.error(f"Transcript file not found: {transcript_file}")
                return context

            logger.debug(f"Reading transcript file: {transcript_file}")
            with open(transcript_file, encoding="utf-8") as tf:
                for line in tf:
                    tr = line.strip()
                    if tr:
                        transcripts.append(tr)

        # Remove duplicates
        transcripts = list(set(transcripts))

        # If no transcripts specified, skip filtering
        if not transcripts:
            logger.debug("No transcripts specified for filtering")
            return context

        logger.info(f"Filtering for {len(transcripts)} transcript(s)")

        # Determine input VCF
        input_vcf = (
            getattr(context, "split_annotations_vcf", None)
            or getattr(context, "filtered_vcf", None)
            or getattr(context, "extracted_vcf", None)
            or context.data
        )

        # Create transcript filter expression
        or_clauses = [f"(EFF[*].TRID = '{tr}')" for tr in transcripts]
        transcript_filter_expr = " | ".join(or_clauses)

        # Generate output filename
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.transcript_filtered.vcf.gz"
        )

        logger.debug(f"Applying transcript filter: {transcript_filter_expr}")

        # Apply transcript filter using SnpSift
        apply_snpsift_filter(
            variant_file=str(input_vcf),
            filter_string=transcript_filter_expr,
            output_file=str(output_vcf),
            cfg=context.config,
        )

        # Update context with filtered VCF
        context.transcript_filtered_vcf = output_vcf  # type: ignore[attr-defined]
        context.data = output_vcf

        logger.info(f"Transcript filtering completed: {output_vcf}")
        return context


class SnpSiftFilterStage(Stage):
    """Apply SnpSift filter expressions."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "snpsift_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Apply SnpSift filter expressions"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Must run after configuration and at least one extraction stage
        # We check for variant extraction explicitly in _process
        return {"configuration_loading"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Run after either extraction stage
        return {"variant_extraction", "parallel_variant_extraction"}

    def _split_before_filter(self, context: PipelineContext) -> bool:
        """Check if we should split before filtering."""
        # Check if snpeff_splitting_mode is set to "before_filters"
        splitting_mode = context.config.get("snpeff_splitting_mode", None)
        return bool(splitting_mode == "before_filters")

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply SnpSift filters."""
        # Check both "filter" and "filters" keys (different parts use different names)
        filter_expr = context.config.get("filter") or context.config.get("filters")

        if not filter_expr or context.config.get("late_filtering"):
            logger.debug("No SnpSift filter or using late filtering")
            return context

        # Get input VCF from context (set by variant extraction stage)
        input_vcf = getattr(context, "extracted_vcf", None) or context.data
        if not input_vcf:
            logger.error("No input VCF available for filtering")
            return context

        # Ensure the input VCF exists
        if not Path(input_vcf).exists():
            logger.error(f"Input VCF file does not exist: {input_vcf}")
            return context

        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.filtered.vcf.gz"
        )

        logger.info(f"Applying SnpSift filter: {filter_expr}")
        filter_config = {"threads": context.config.get("threads", 1)}
        apply_snpsift_filter(str(input_vcf), filter_expr, filter_config, str(output_vcf))

        context.filtered_vcf = output_vcf
        # Only update context.data if we haven't extracted to TSV yet
        if not hasattr(context, "extracted_tsv"):
            context.data = output_vcf

        return context


class FieldExtractionStage(Stage):
    """Extract specified fields from VCF to TSV format."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "field_extraction"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Extract fields from VCF to TSV format"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Must have gene bed and config, and variant extraction
        deps = {"gene_bed_creation", "configuration_loading", "variant_extraction"}
        return deps

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # Prefer to run after filtering if it exists
        return {"snpsift_filtering", "transcript_filtering"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract fields to TSV format."""
        # Use the most recent VCF in the processing chain
        input_vcf = (
            getattr(context, "transcript_filtered_vcf", None)
            or getattr(context, "filtered_vcf", None)
            or getattr(context, "split_annotations_vcf", None)
            or getattr(context, "extracted_vcf", None)
            or context.data
        )

        logger.debug(f"Using VCF for field extraction: {input_vcf}")

        # Get fields_to_extract as string from config (not "extract")
        fields_str = context.config.get("fields_to_extract", "")

        # Handle extra sample fields if requested - append to the fields string
        if context.config.get("append_extra_sample_fields", False) and context.config.get(
            "extra_sample_fields", []
        ):
            fields_str = ensure_fields_in_extract(fields_str, context.config["extra_sample_fields"])
            logger.debug(f"Updated fields with extra sample fields: {fields_str}")

        # Convert string to list for processing
        fields = fields_str.split() if fields_str else []

        # Debug logging
        logger.debug(f"FieldExtractionStage - config keys: {list(context.config.keys())[:10]}...")
        logger.debug(f"FieldExtractionStage - 'extract' in config: {'extract' in context.config}")
        logger.debug(
            "FieldExtractionStage - 'fields_to_extract' in config: "
            f"{'fields_to_extract' in context.config}"
        )
        logger.debug(f"FieldExtractionStage - fields value: {fields}")

        if not fields:
            raise ValueError("No fields specified for extraction")

        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.tsv"
        )

        # Handle gzip compression if requested (default: True)
        if context.config.get("gzip_intermediates", True):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Extracting {len(fields)} fields from VCF to TSV using bcftools query")

        # Prepare config for extract_fields_bcftools
        extract_config = {
            "debug_level": context.config.get("log_level", "INFO"),
        }

        # Convert fields list to space-separated string
        fields_str = " ".join(fields)

        # Use bcftools query (19x faster than SnpSift)
        extract_fields_bcftools(
            variant_file=str(input_vcf),
            fields=fields_str,
            cfg=extract_config,
            output_file=str(output_tsv),
            vcf_samples=context.vcf_samples,  # For per-sample column naming
        )

        context.extracted_tsv = output_tsv
        context.data = output_tsv

        return context

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the extracted TSV file."""
        if context.extracted_tsv:
            return [context.extracted_tsv]
        return []


class GenotypeReplacementStage(Stage):
    """Replace genotype values with sample IDs."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "genotype_replacement"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Replace genotype values with sample IDs"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"field_extraction", "sample_config_loading"}

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the context paths
        so that subsequent stages can find the genotype_replaced_tsv file.
        """
        # Reconstruct the expected output path
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.genotype_replaced.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates", True):
            output_tsv = Path(str(output_tsv) + ".gz")

        # Check if the file exists and restore context
        if output_tsv.exists():
            context.genotype_replaced_tsv = output_tsv
            context.data = output_tsv
            logger.info(
                f"Restored genotype_replaced_tsv path to {output_tsv} after checkpoint skip"
            )
        else:
            logger.warning(f"Expected genotype_replaced_tsv file not found: {output_tsv}")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """No-op: genotype replacement eliminated in Phase 11.

        GT reconstruction is deferred to output time (TSV/Excel stages).
        Raw per-sample GT columns from bcftools flow directly to analysis.
        """
        logger.info(
            "Genotype replacement skipped — raw per-sample GT columns flow to analysis (Phase 11)"
        )
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        if context.data:
            return [context.data]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            return [context.genotype_replaced_tsv]
        return []


class PhenotypeIntegrationStage(Stage):
    """Integrate phenotype data into the variant table."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_integration"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Add phenotype data to variant table"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Only declare the hard requirement - field extraction must be complete
        return {"field_extraction"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        # These stages should run before if they exist in the pipeline
        return {"genotype_replacement", "phenotype_loading"}

    def _has_genotype_replacement(self, context: PipelineContext) -> bool:
        """Check if genotype replacement is enabled."""
        # Check if genotype replacement is not disabled and GT field is present
        if context.config.get("no_replacement", False):
            return False

        # Check if genotype_replaced_tsv exists in context
        return (
            hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv is not None
        )

    def _has_phenotype_data(self, context: PipelineContext) -> bool:
        """Check if phenotype data is available."""
        # Check if phenotypes were loaded
        if not hasattr(context, "phenotype_data") or context.phenotype_data is None:
            return False

        # Check if phenotypes dict is populated
        return bool(context.phenotype_data)

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to restore the context paths
        so that subsequent stages can find the phenotypes_added_tsv file.
        """
        # Reconstruct the expected output path
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.phenotypes_added.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        # Check if the file exists and restore context
        if output_tsv.exists():
            context.phenotypes_added_tsv = output_tsv
            context.data = output_tsv
            logger.info(f"Restored phenotypes_added_tsv path to {output_tsv} after checkpoint skip")
        else:
            logger.warning(f"Expected phenotypes_added_tsv file not found: {output_tsv}")

        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Add phenotype data to the table."""
        if not context.phenotype_data:
            logger.debug("No phenotype data to integrate")
            return context

        # Determine input file - use genotype_replaced if available, otherwise extracted
        if hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            input_tsv = context.genotype_replaced_tsv
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            input_tsv = context.extracted_tsv
        else:
            input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.phenotypes_added.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Adding phenotype data for {len(context.vcf_samples)} samples")

        # Read input TSV
        compression = "gzip" if str(input_tsv).endswith(".gz") else None
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=compression)

        # Phase 11: Check for per-sample columns (bcftools output) vs packed GT column
        # Per-sample columns are sample IDs directly (from context.vcf_samples)
        has_sample_columns = context.vcf_samples and all(
            sample in df.columns for sample in context.vcf_samples
        )

        if has_sample_columns:
            # Phase 11: Use per-sample columns directly
            logger.debug("Using per-sample GT columns for phenotype extraction")
            df["Phenotypes"] = df.apply(
                lambda row: extract_phenotypes_from_sample_columns(
                    row, context.vcf_samples, context.phenotype_data
                ),
                axis=1,
            )
        elif "GT" in df.columns:
            # Backwards compatibility: Use packed GT column (from genotype replacement)
            logger.debug("Using packed GT column for phenotype extraction (legacy mode)")
            df["Phenotypes"] = df["GT"].apply(
                lambda gt_val: extract_phenotypes_for_gt_row(gt_val, context.phenotype_data)
            )
        else:
            logger.warning("No GT column or per-sample columns found, cannot add phenotype data")
            return context

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=compression)

        context.phenotypes_added_tsv = output_tsv
        context.data = output_tsv
        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        if hasattr(context, "genotype_replaced_tsv") and context.genotype_replaced_tsv:
            return [context.genotype_replaced_tsv]
        elif hasattr(context, "extracted_tsv") and context.extracted_tsv:
            return [context.extracted_tsv]
        return []

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "phenotypes_added_tsv") and context.phenotypes_added_tsv:
            return [context.phenotypes_added_tsv]
        return []


class ExtraColumnRemovalStage(Stage):
    """Remove extra columns from the TSV file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "extra_column_removal"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Remove specified columns from TSV"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"phenotype_integration", "genotype_replacement", "field_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Remove extra columns if specified."""
        columns_to_remove = context.config.get("extra_columns_to_remove", [])

        if not columns_to_remove:
            logger.debug("No columns to remove")
            return context

        input_tsv = context.data
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.columns_removed.tsv"
        )

        # Handle gzip compression
        if context.config.get("gzip_intermediates"):
            output_tsv = Path(str(output_tsv) + ".gz")

        logger.info(f"Removing {len(columns_to_remove)} columns")

        # Read, remove columns, and write
        compression = "gzip" if str(input_tsv).endswith(".gz") else None
        df = pd.read_csv(input_tsv, sep="\t", dtype=str, compression=compression)

        # Remove columns that exist
        existing_cols = [col for col in columns_to_remove if col in df.columns]
        if existing_cols:
            df = df.drop(columns=existing_cols)
            logger.info(f"Removed columns: {existing_cols}")

        # Write output
        compression = "gzip" if str(output_tsv).endswith(".gz") else None
        df.to_csv(output_tsv, sep="\t", index=False, compression=compression)

        context.extra_columns_removed_tsv = output_tsv
        context.data = output_tsv
        return context


class StreamingDataProcessingStage(Stage):
    """Single-pass streaming processing for memory efficiency."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "streaming_processing"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Stream process large files in a single pass"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        # Can replace the separate genotype/phenotype/column stages
        return {"field_extraction"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return False  # Streaming is inherently sequential

    @property
    def memory_efficient(self) -> bool:
        """Return whether this stage processes data in a memory-efficient manner."""
        return True  # Processes file line-by-line without loading entire content

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process file in streaming fashion for memory efficiency."""
        # Check if we should use streaming
        if not context.config.get("use_streaming_processing"):
            logger.debug("Streaming processing not enabled")
            return context

        input_file = context.data
        output_file = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.processed.tsv"
        )

        # Handle compression
        if context.config.get("gzip_intermediates"):
            output_file = Path(str(output_file) + ".gz")

        logger.info("Processing file in streaming mode for memory efficiency")

        # Get processing parameters
        replace_genotypes = context.config.get("replace_genotypes", False)
        samples = context.vcf_samples if replace_genotypes else []
        phenotypes = context.phenotype_data if context.phenotype_data else {}
        columns_to_remove = context.config.get("extra_columns_to_remove", [])

        # Stream process the file
        self._stream_process_file(
            input_file=input_file,
            output_file=output_file,
            samples=samples,
            phenotypes=phenotypes,
            columns_to_remove=columns_to_remove,
            missing_string=context.config.get("missing_string", "./."),
        )

        context.data = output_file
        return context

    def _stream_process_file(
        self,
        input_file: Path,
        output_file: Path,
        samples: list[str],
        phenotypes: dict[str, str],
        columns_to_remove: list[str],
        missing_string: str,
    ) -> None:
        """Stream process file line by line."""

        # Validate input file
        input_file = validate_file_exists(input_file, self.name)

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Determine compression
        in_compressed = str(input_file).endswith(".gz")
        out_compressed = str(output_file).endswith(".gz")

        try:
            # Open files with proper error handling
            try:
                in_fh_obj = gzip.open(input_file, "rt") if in_compressed else open(input_file)  # noqa: SIM115 - used in `with` below, try/except needed for error mapping
            except OSError as e:
                raise FileFormatError(str(input_file), "readable TSV file", self.name) from e

            try:
                out_fh_obj = (
                    gzip.open(output_file, "wt") if out_compressed else open(output_file, "w")  # noqa: SIM115 - used in `with` below
                )
            except OSError as e:
                in_fh_obj.close()
                raise PermissionError(f"Cannot write to output file: {output_file}") from e

            with in_fh_obj as in_fh, out_fh_obj as out_fh:
                # Process header
                header = in_fh.readline().strip()
                if not header:
                    raise FileFormatError(str(input_file), "TSV file with header", self.name)

                columns = header.split("\t")

                # Find columns to process
                gt_columns = [i for i, col in enumerate(columns) if col.endswith(".GT")]
                columns_to_remove_idx = [
                    i for i, col in enumerate(columns) if col in columns_to_remove
                ]

                # Update header
                new_columns = []
                for i, col in enumerate(columns):
                    if i in columns_to_remove_idx:
                        continue
                    if samples and i < len(columns) and columns[i].endswith(".GT"):
                        # Replace with sample name
                        sample_idx = gt_columns.index(i)
                        if sample_idx < len(samples):
                            new_columns.append(samples[sample_idx])
                        else:
                            new_columns.append(col)
                    else:
                        new_columns.append(col)

                # Add phenotype column if needed
                if phenotypes:
                    new_columns.append("Phenotypes")

                # Write header
                out_fh.write("\t".join(new_columns) + "\n")

                # Process data lines
                for line in in_fh:
                    fields = line.strip().split("\t")

                    # Build new row
                    new_fields = []
                    for i, field in enumerate(fields):
                        if i in columns_to_remove_idx:
                            continue

                        if samples and i in gt_columns:
                            # Replace genotype encoding
                            if field in ("0/1", "1/0", "0|1", "1|0"):
                                new_fields.append("het")
                            elif field in ["1/1", "1|1"]:
                                new_fields.append("hom")
                            elif field in ["0/0", "0|0"]:
                                new_fields.append("ref")
                            elif field == missing_string:
                                new_fields.append("missing")
                            else:
                                new_fields.append(field)
                        else:
                            new_fields.append(field)

                    # Add phenotypes if needed
                    if phenotypes:
                        pheno_values = []
                        for sample in samples:
                            pheno_values.append(phenotypes.get(sample, ""))
                        new_fields.append(";".join(pheno_values))

                    # Write line
                    out_fh.write("\t".join(new_fields) + "\n")

        except Exception as e:
            logger.error(f"Error during streaming processing: {e}")
            # Clean up partial output file
            if output_file.exists():
                with contextlib.suppress(Exception):
                    output_file.unlink()
            raise

        logger.info(f"Streaming processing complete: {output_file}")


class ParallelCompleteProcessingStage(Stage):
    """Run complete processing pipeline in parallel for each BED chunk.

    This stage replicates the old pipeline's behavior where each chunk
    is processed completely (extraction -> filtering -> field extraction)
    in parallel, then the TSV results are merged.
    """

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "parallel_complete_processing"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Process variants in parallel with complete pipeline per chunk"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"gene_bed_creation", "configuration_loading"}

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        # This stage manages its own parallelism internally
        return False

    @property
    def estimated_runtime(self) -> float:
        """Return the estimated runtime in seconds."""
        return 60.0  # Complete pipeline takes longer

    def _validate_chunk_tsv(self, chunk_tsv: Path, fallback_on_error: bool = True) -> bool:
        """Validate that a chunk TSV file is complete and valid.

        Parameters
        ----------
        chunk_tsv : Path
            Path to the TSV chunk file to validate
        fallback_on_error : bool
            If True, validation errors result in False (file will be reprocessed)
            If False, validation errors raise exceptions

        Returns
        -------
        bool
            True if file is valid, False if invalid or missing
        """
        if not chunk_tsv.exists():
            logger.debug(f"Chunk TSV {chunk_tsv} does not exist")
            return False

        try:
            # Size check (must be > 0)
            stat = chunk_tsv.stat()
            if stat.st_size == 0:
                logger.debug(f"Chunk TSV {chunk_tsv} is empty, invalid")
                return False

            # Basic content validation for TSV files
            try:
                import gzip

                open_func = gzip.open if chunk_tsv.name.endswith(".gz") else open
                mode = "rt" if chunk_tsv.name.endswith(".gz") else "r"

                with open_func(chunk_tsv, mode) as f:
                    first_line = f.readline().strip()
                    # Check if it looks like a TSV header
                    if "\t" not in first_line:
                        logger.debug(f"Chunk TSV {chunk_tsv} has invalid header format")
                        return False
                    # Try to read one more line to ensure file is not just header
                    second_line = f.readline()
                    if not second_line:
                        logger.debug(f"Chunk TSV {chunk_tsv} only has header, no data")
                        return False
            except Exception as e:
                logger.debug(f"Chunk TSV {chunk_tsv} failed content validation: {e}")
                if not fallback_on_error:
                    raise
                return False

            logger.debug(f"Chunk TSV {chunk_tsv} is valid (size: {stat.st_size} bytes)")
            return True

        except Exception as e:
            logger.warning(f"Error validating chunk TSV {chunk_tsv}: {e}")
            if not fallback_on_error:
                raise
            return False

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        When this stage is skipped, we need to validate chunk-level outputs
        and mark the constituent stages as complete so that dependent stages can proceed.
        """
        # Use the workspace base_name to reconstruct the expected output file path
        base_name = context.workspace.base_name
        intermediate_dir = context.workspace.intermediate_dir

        # Set the expected merged TSV path
        merged_tsv = intermediate_dir / f"{base_name}.extracted.tsv"

        # Check if gzipped version exists (which is more likely with compression)
        if (intermediate_dir / f"{base_name}.extracted.tsv.gz").exists():
            merged_tsv = intermediate_dir / f"{base_name}.extracted.tsv.gz"

        # Validate the merged TSV exists and is valid
        if merged_tsv.exists() and self._validate_chunk_tsv(merged_tsv):
            logger.info(f"Restored valid merged TSV: {merged_tsv}")
        else:
            logger.warning(f"Expected merged TSV not found or invalid: {merged_tsv}")

        # Restore context state that would have been set
        context.extracted_tsv = merged_tsv
        context.data = merged_tsv

        # Mark all the stages this stage would have completed
        context.mark_complete("variant_extraction")
        context.mark_complete("parallel_variant_extraction")
        context.mark_complete("snpsift_filtering")
        context.mark_complete("field_extraction")

        logger.info(
            f"Restored data path to {merged_tsv} and marked constituent stages as complete "
            f"after checkpoint skip"
        )
        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process variants in parallel with complete pipeline per chunk."""
        threads = context.config.get("threads", 1)

        if threads <= 1:
            logger.info("Thread count <= 1, falling back to sequential stages")
            # Return context unchanged - let normal stages handle it
            return context

        # Split BED file into chunks
        if context.gene_bed_file is None:
            raise ValueError("gene_bed_file is required for parallel complete processing")
        split_start = self._start_subtask("bed_splitting")

        bed_chunks = self._split_bed_file(context.gene_bed_file, threads)
        self._end_subtask("bed_splitting", split_start)
        logger.info(f"Split BED file into {len(bed_chunks)} chunks for parallel processing")

        # Process chunks in parallel (complete pipeline per chunk)
        process_start = self._start_subtask("parallel_chunk_processing")
        chunk_tsvs = self._process_chunks_parallel(context, bed_chunks)
        self._end_subtask("parallel_chunk_processing", process_start)

        # Merge TSV results
        merge_start = self._start_subtask("tsv_merging")
        merged_tsv = self._merge_tsv_outputs(context, chunk_tsvs)
        self._end_subtask("tsv_merging", merge_start)

        # Update context
        context.extracted_tsv = merged_tsv
        context.data = merged_tsv

        # Mark stages as complete so they don't run again
        context.mark_complete("variant_extraction")
        context.mark_complete("parallel_variant_extraction")
        context.mark_complete("snpsift_filtering")
        context.mark_complete("field_extraction")

        # Mark chunked processing as complete so analysis stages can run
        context.config["chunked_processing_complete"] = True

        # Cleanup chunks
        cleanup_start = self._start_subtask("cleanup")
        self._cleanup_chunks(bed_chunks, chunk_tsvs, context)
        self._end_subtask("cleanup", cleanup_start)

        return context

    def get_input_files(self, context: PipelineContext) -> list[Path]:
        """Return input files for checkpoint tracking."""
        input_files = []
        if hasattr(context, "gene_bed_file") and context.gene_bed_file:
            input_files.append(context.gene_bed_file)
        if hasattr(context, "vcf_file") and context.vcf_file:
            input_files.append(Path(context.vcf_file))
        return input_files

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return output files for checkpoint tracking."""
        if hasattr(context, "extracted_tsv") and context.extracted_tsv:
            return [context.extracted_tsv]
        return []

    def _split_bed_file(self, bed_file: Path, n_chunks: int) -> list[Path]:
        """Split BED file into roughly equal chunks."""
        chunk_dir = bed_file.parent / "chunks"
        chunk_dir.mkdir(exist_ok=True)

        chunks = split_bed_file(str(bed_file), n_chunks, str(chunk_dir))
        return [Path(chunk) for chunk in chunks]

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
        import time

        chunk_base = f"{base_name}.chunk_{chunk_index}"

        # Step 1: Extract variants
        extract_start = time.time()
        chunk_vcf = intermediate_dir / f"{chunk_base}.variants.vcf.gz"
        extract_variants(
            vcf_file=vcf_file,
            bed_file=str(chunk_bed),
            cfg={
                "threads": config.get("threads_per_chunk", 2),
                "bcftools_prefilter": config.get("bcftools_prefilter"),
            },
            output_file=str(chunk_vcf),
        )
        extract_time = time.time() - extract_start

        # Step 2: Apply SnpSift filter (if not late filtering)
        filter_start = time.time()
        if config.get("late_filtering", False):
            filtered_vcf = chunk_vcf
            filter_time = 0.0
        else:
            filtered_vcf = intermediate_dir / f"{chunk_base}.filtered.vcf.gz"
            filter_expr = config.get("filter") or config.get("filters")
            if filter_expr:
                apply_snpsift_filter(
                    str(chunk_vcf),
                    filter_expr,
                    {"threads": config.get("threads_per_chunk", 2)},
                    str(filtered_vcf),
                )
            else:
                filtered_vcf = chunk_vcf
            filter_time = time.time() - filter_start

        # Step 3: Extract fields to TSV
        field_extract_start = time.time()
        # Respect gzip_intermediates setting for compression
        use_compression = config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"
        chunk_tsv = intermediate_dir / f"{chunk_base}{tsv_suffix}"

        # Get fields_to_extract as string from config (not "extract")
        fields_str = config.get("fields_to_extract", "")

        # Handle extra sample fields if requested - append to the fields string
        if config.get("append_extra_sample_fields", False) and config.get(
            "extra_sample_fields", []
        ):
            fields_str = ensure_fields_in_extract(fields_str, config["extra_sample_fields"])
            logger.debug(f"Updated fields with extra sample fields: {fields_str}")

        # Convert string to list for processing
        fields = fields_str.split() if fields_str else []
        if not fields:
            raise ValueError("No fields specified for extraction")

        extract_fields_bcftools(
            variant_file=str(filtered_vcf),
            fields=" ".join(fields),
            cfg={"extract_fields_separator": config.get("extract_fields_separator", ",")},
            output_file=str(chunk_tsv),
            vcf_samples=config.get("vcf_samples"),
        )
        field_extract_time = time.time() - field_extract_start

        # Store subtask times if dict provided
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
        vcf_file = context.config["vcf_file"]
        threads = context.config.get("threads", 1)
        base_name = context.workspace.base_name
        intermediate_dir = context.workspace.intermediate_dir

        # Handle empty bed_chunks case
        if not bed_chunks:
            logger.warning("No BED chunks to process - skipping parallel processing")
            return []

        # Phase 1: Check for existing valid chunk TSVs
        chunks_to_process = []
        existing_chunk_tsvs = []

        use_compression = context.config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"

        for i, chunk_bed in enumerate(bed_chunks):
            chunk_base = f"{base_name}.chunk_{i}"
            expected_tsv = intermediate_dir / f"{chunk_base}{tsv_suffix}"

            if self._validate_chunk_tsv(expected_tsv):
                logger.info(f"Reusing existing chunk TSV {i}: {expected_tsv.name}")
                existing_chunk_tsvs.append(expected_tsv)
            else:
                logger.info(f"Processing missing/invalid chunk TSV {i}")
                chunks_to_process.append((i, chunk_bed, expected_tsv))

        # Phase 2: Process only missing chunks
        new_chunk_tsvs = []
        if chunks_to_process:
            # Each worker gets limited threads
            threads_per_worker = max(2, threads // len(chunks_to_process))

            # Prepare config for workers - ensure field extraction works correctly
            worker_config = {
                "threads_per_chunk": threads_per_worker,
                "bcftools_prefilter": context.config.get("bcftools_prefilter"),
                "late_filtering": context.config.get("late_filtering", False),
                "filter": context.config.get("filter"),
                "filters": context.config.get("filters"),
                "fields_to_extract": context.config.get("fields_to_extract", ""),
                "append_extra_sample_fields": context.config.get(
                    "append_extra_sample_fields", False
                ),
                "extra_sample_fields": context.config.get("extra_sample_fields", []),
                "extract_fields_separator": context.config.get("extract_fields_separator", ","),
                "gzip_intermediates": use_compression,
                "vcf_samples": context.vcf_samples,
            }

            # Use a manager to share the subtask times dict across processes
            from multiprocessing import Manager

            manager = Manager()
            shared_subtask_times = manager.dict()

            logger.info(
                f"Processing {len(chunks_to_process)} missing chunks out of {len(bed_chunks)} total"
            )

            with ProcessPoolExecutor(max_workers=len(chunks_to_process)) as executor:
                # Submit jobs for missing chunks only
                future_to_chunk = {}
                for i, chunk_bed, expected_tsv in chunks_to_process:
                    future = executor.submit(
                        self._process_single_chunk,
                        i,
                        chunk_bed,
                        vcf_file,
                        base_name,
                        intermediate_dir,
                        worker_config,
                        shared_subtask_times,
                    )
                    future_to_chunk[future] = (i, chunk_bed, expected_tsv)

                # Collect results
                for future in as_completed(future_to_chunk):
                    chunk_index, chunk_bed, expected_tsv = future_to_chunk[future]
                    try:
                        chunk_tsv = future.result()
                        new_chunk_tsvs.append(chunk_tsv)
                        logger.debug(f"Completed chunk {chunk_index}: {chunk_bed.name}")
                    except Exception as e:
                        logger.error(f"Failed to process chunk {chunk_bed}: {e}")
                        raise

            # Aggregate subtask times from new chunks
            if shared_subtask_times:
                self._aggregate_subtask_times(shared_subtask_times)
        else:
            logger.info("All chunk TSVs already exist and are valid - no processing needed")

        # Phase 3: Combine existing and new outputs, maintaining order
        all_chunk_tsvs = []
        for i in range(len(bed_chunks)):
            chunk_base = f"{base_name}.chunk_{i}"
            expected_tsv = intermediate_dir / f"{chunk_base}{tsv_suffix}"
            all_chunk_tsvs.append(expected_tsv)

        return all_chunk_tsvs

    def _aggregate_subtask_times(self, shared_subtask_times: Any) -> None:
        """Aggregate subtask times from parallel processing."""
        # Group by subtask type
        extraction_times = []
        filtering_times = []
        field_extraction_times = []

        for key, value in shared_subtask_times.items():
            if "_extraction" in key and "_field_" not in key:
                extraction_times.append(value)
            elif "_filtering" in key:
                filtering_times.append(value)
            elif "_field_extraction" in key:
                field_extraction_times.append(value)

        # Record aggregated times
        if extraction_times:
            self._subtask_times["variant_extraction"] = sum(extraction_times) / len(
                extraction_times
            )
        if filtering_times:
            self._subtask_times["snpsift_filtering"] = sum(filtering_times) / len(filtering_times)
        if field_extraction_times:
            self._subtask_times["field_extraction"] = sum(field_extraction_times) / len(
                field_extraction_times
            )

    def _merge_tsv_outputs(self, context: PipelineContext, chunk_tsvs: list[Path]) -> Path:
        """Merge TSV chunks into a single file with optimized compression handling."""
        output_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.tsv"
        )

        # Always use gzip for intermediate files to save space
        use_compression = context.config.get("gzip_intermediates", True)
        if use_compression:
            output_tsv = Path(str(output_tsv) + ".gz")

        if len(chunk_tsvs) == 1:
            # Just move the single file
            shutil.move(str(chunk_tsvs[0]), str(output_tsv))
            logger.info(f"Moved single TSV chunk to {output_tsv}")
        else:
            # Use optimized merging based on compression
            if use_compression and all(str(chunk).endswith(".gz") for chunk in chunk_tsvs):
                self._merge_gzipped_tsvs_optimized(chunk_tsvs, output_tsv)
            else:
                self._merge_tsvs_traditional(chunk_tsvs, output_tsv)

            logger.info(f"Merged {len(chunk_tsvs)} TSV chunks into {output_tsv}")

        return output_tsv

    def _merge_gzipped_tsvs_optimized(self, chunk_tsvs: list[Path], output_tsv: Path) -> None:
        """Optimized merging of gzipped TSV files using direct concatenation where possible."""
        # Sort chunks to ensure consistent order
        sorted_chunks = sorted(chunk_tsvs)

        # Try to use fast concatenation for identical headers
        if self._can_use_fast_concatenation(sorted_chunks):
            logger.info("Using optimized gzipped file concatenation")
            self._fast_concatenate_gzipped_tsvs(sorted_chunks, output_tsv)
        else:
            logger.info("Using traditional merging due to header differences")
            self._merge_tsvs_traditional(sorted_chunks, output_tsv)

    def _can_use_fast_concatenation(self, chunk_tsvs: list[Path]) -> bool:
        """Check if all chunks have identical headers for fast concatenation."""
        if len(chunk_tsvs) <= 1:
            return True

        try:
            # Read first header
            with gzip.open(chunk_tsvs[0], "rt") as f:
                first_header = f.readline().strip()

            # Check all other headers
            for chunk in chunk_tsvs[1:]:
                with gzip.open(chunk, "rt") as f:
                    header = f.readline().strip()
                    if header != first_header:
                        return False
            return True
        except Exception as e:
            logger.debug(f"Error checking headers for fast concatenation: {e}")
            return False

    def _fast_concatenate_gzipped_tsvs(self, chunk_tsvs: list[Path], output_tsv: Path) -> None:
        """Fast concatenation of gzipped TSVs using streaming compression."""

        # Use fastest compression level for intermediate files
        compression_level = 1  # Fastest compression

        with gzip.open(output_tsv, "wt", compresslevel=compression_level) as out_fh:
            first_file = True

            for chunk_tsv in chunk_tsvs:
                with gzip.open(chunk_tsv, "rt") as in_fh:
                    if first_file:
                        # Copy entire first file including header
                        shutil.copyfileobj(in_fh, out_fh)
                        first_file = False
                    else:
                        # Skip header line and copy rest
                        next(in_fh, None)  # Skip header
                        shutil.copyfileobj(in_fh, out_fh)

    def _merge_tsvs_traditional(self, chunk_tsvs: list[Path], output_tsv: Path) -> None:
        """Traditional TSV merging with compression handling."""

        # Open output file
        with (
            gzip.open(output_tsv, "wt", compresslevel=1)
            if str(output_tsv).endswith(".gz")
            else open(output_tsv, "w")
        ) as out_fh:
            first_file = True
            for chunk_tsv in sorted(chunk_tsvs):  # Sort to ensure consistent order
                # Open chunk file
                with (
                    gzip.open(chunk_tsv, "rt")
                    if str(chunk_tsv).endswith(".gz")
                    else open(chunk_tsv)
                ) as in_fh:
                    if first_file:
                        # Copy entire first file including header
                        shutil.copyfileobj(in_fh, out_fh)
                        first_file = False
                    else:
                        # Skip header line
                        next(in_fh, None)
                        # Copy rest of file
                        shutil.copyfileobj(in_fh, out_fh)

    def _cleanup_chunks(
        self, bed_chunks: list[Path], tsv_chunks: list[Path], context: PipelineContext
    ) -> None:
        """Clean up temporary chunk files."""
        if not context.config.get("keep_intermediates", False):
            for chunk in bed_chunks + tsv_chunks:
                if chunk.exists():
                    chunk.unlink()

            # Remove chunks directory if empty
            if bed_chunks and bed_chunks[0].parent.name == "chunks":
                chunk_dir = bed_chunks[0].parent
                if chunk_dir.exists() and not any(chunk_dir.iterdir()):
                    chunk_dir.rmdir()


class DataSortingStage(Stage):
    """Sort extracted TSV data by gene column for efficient chunked processing."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "data_sorting"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Sort TSV data by gene column for efficient processing"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"field_extraction"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Sort the extracted TSV by gene column if gene column is present."""
        if not context.extracted_tsv:
            logger.warning("No extracted TSV found for sorting")
            return context

        # Check if sorting is needed and if gene column exists
        gene_column = context.config.get("gene_column_name", "GENE")

        # First check if gene column exists in the TSV
        has_gene_column = False
        try:
            with smart_open(str(context.extracted_tsv), "r") as f:
                header = f.readline().strip()
                if header:
                    columns = header.split("\t")
                    has_gene_column = gene_column in columns
        except Exception as e:
            logger.warning(f"Could not read header from {context.extracted_tsv}: {e}")
            return context

        if not has_gene_column:
            logger.info(f"Gene column '{gene_column}' not found in TSV, skipping sorting")
            return context

        # Create sorted output file
        sorted_tsv = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.extracted.sorted.tsv"
        )

        # Handle gzip compression if requested
        if context.config.get("gzip_intermediates"):
            sorted_tsv = Path(str(sorted_tsv) + ".gz")

        logger.info(f"Sorting TSV by gene column '{gene_column}' for efficient processing")

        # Use the sort function from old pipeline
        self._sort_tsv_by_gene(
            input_file=str(context.extracted_tsv),
            output_file=str(sorted_tsv),
            gene_column=gene_column,
            temp_dir=str(context.workspace.intermediate_dir),
            memory_limit=context.config.get("sort_memory_limit", "2G"),
            parallel=context.config.get("sort_parallel", 4),
        )

        # Clean up unsorted file if not keeping intermediates
        if not context.config.get("keep_intermediates", False):
            with contextlib.suppress(OSError):
                os.remove(str(context.extracted_tsv))

        # Update context with sorted file
        context.extracted_tsv = sorted_tsv
        context.data = sorted_tsv

        return context

    def _sort_tsv_by_gene(
        self,
        input_file: str,
        output_file: str,
        gene_column: str = "GENE",
        temp_dir: str | None = None,
        memory_limit: str = "2G",
        parallel: int = 4,
    ) -> str:
        """
        Sort a TSV file by gene column using external sort command for memory efficiency.

        This is adapted from the original pipeline's sort_tsv_by_gene function.
        """
        logger.info(f"Sorting TSV file by gene column: {input_file} -> {output_file}")
        logger.info(f"Using memory limit: {memory_limit}, parallel threads: {parallel}")

        # First, find the gene column index
        with smart_open(input_file, "r") as f:
            header = f.readline().strip()
            if not header:
                raise ValueError(f"Input file '{input_file}' is empty or has no header")
            columns = header.split("\t")

        try:
            gene_col_idx = columns.index(gene_column) + 1  # sort uses 1-based indexing
        except ValueError as e:
            raise ValueError(f"Gene column '{gene_column}' not found in TSV file") from e

        logger.info(f"Found gene column '{gene_column}' at position {gene_col_idx}")

        # Resolve full paths for shell utilities to avoid PATH issues in
        # subprocess shells (e.g., WSL, conda envs, non-interactive bash).
        import shlex

        gzip_bin = shutil.which("gzip") or "gzip"
        sort_bin = shutil.which("sort") or "sort"
        tail_bin = shutil.which("tail") or "tail"

        # Build sort arguments with memory efficiency options
        sort_args = [
            f"-k{gene_col_idx},{gene_col_idx}",  # Sort by gene column
            f"--buffer-size={memory_limit}",  # Memory limit
            f"--parallel={parallel}",  # Parallel threads
            "--stable",  # Stable sort for consistent results
            f"--compress-program={shlex.quote(gzip_bin)}",  # Use gzip for temp files
        ]

        if temp_dir:
            sort_args.extend(["-T", temp_dir])

        sort_cmd = " ".join(sort_args)

        safe_input = shlex.quote(input_file)
        safe_output = shlex.quote(output_file)

        # Build the complete command based on compression
        if input_file.endswith(".gz"):
            if output_file.endswith(".gz"):
                # Both gzipped
                cmd = (
                    f"{gzip_bin} -cd {safe_input} | "
                    f'{{ IFS= read -r header; echo "$header"; '
                    f"{sort_bin} {sort_cmd} -t$'\\t'; }} | "
                    f"{gzip_bin} -c > {safe_output}"
                )
            else:
                # Input gzipped, output not
                cmd = (
                    f"{gzip_bin} -cd {safe_input} | "
                    f"{{ IFS= read -r header; echo \"$header\"; {sort_bin} {sort_cmd} -t$'\\t'; }} "
                    f"> {safe_output}"
                )
        else:
            if output_file.endswith(".gz"):
                # Input not gzipped, output gzipped
                cmd = (
                    f'{{ IFS= read -r header < {safe_input}; echo "$header"; '
                    f"{tail_bin} -n +2 {safe_input} | {sort_bin} {sort_cmd} -t$'\\t'; }} | "
                    f"{gzip_bin} -c > {safe_output}"
                )
            else:
                # Neither gzipped
                cmd = (
                    f'{{ IFS= read -r header < {safe_input}; echo "$header"; '
                    f"{tail_bin} -n +2 {safe_input} | {sort_bin} {sort_cmd} -t$'\\t'; }} "
                    f"> {safe_output}"
                )

        # Execute the command
        logger.debug(f"Running sort command: {cmd}")

        # Find bash executable
        bash_candidate = shutil.which("bash") or "/bin/bash"
        bash_path: str | None = bash_candidate if os.path.exists(bash_candidate) else None
        if bash_path is None:
            logger.warning(f"bash not found at {bash_candidate}, falling back to shell default")

        # Execute command
        result = subprocess.run(
            cmd, shell=True, capture_output=True, text=True, executable=bash_path
        )

        if result.returncode != 0:
            logger.error(f"Sort command failed: {result.stderr}")
            raise RuntimeError(f"Failed to sort TSV file: {result.stderr}")

        # Verify output file was created
        if not os.path.exists(output_file):
            raise RuntimeError(f"Output file '{output_file}' was not created")

        # Check if output file has content
        with smart_open(output_file, "r") as f:
            first_line = f.readline()
            if not first_line:
                raise RuntimeError(f"Output file '{output_file}' is empty")

        logger.info("Successfully sorted TSV file by gene column")
        return output_file

    def get_output_files(self, context: PipelineContext) -> list[Path]:
        """Return the sorted TSV file."""
        if context.extracted_tsv:
            return [context.extracted_tsv]
        return []

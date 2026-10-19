"""
Processing stages for variant extraction and data transformation.

This module contains stages that handle the core data processing tasks:
- Gene BED file creation
- Variant extraction from VCF files
- Filtering with bcftools and SnpSift
- Field extraction to TSV format
- Phenotype integration
"""

import gzip
import logging
import os
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd  # noqa: F401
from smart_open import smart_open  # noqa: F401

from ..extractor import extract_fields_bcftools, get_expected_extraction_columns
from ..filters import apply_snpsift_filter, extract_variants
from ..gene_bed import get_gene_bed, normalize_genes
from ..pipeline_core import PipelineContext, Stage
from ..pipeline_core.error_handling import (
    FileFormatError,
    ToolNotFoundError,
    graceful_error_handling,
    retry_on_failure,
    validate_file_exists,
)
from ..transcript_filter import filter_vcf_to_transcripts, load_transcript_ids
from ..utils import ensure_fields_in_extract, split_bed_file
from ..vcf_eff_one_per_line import process_vcf_file as split_snpeff_annotations
from .post_processing_stages import (
    DataSortingStage,
    ExtraColumnRemovalStage,
    GenotypeReplacementStage,
    PCAComputationStage,
    PhenotypeIntegrationStage,
)

logger = logging.getLogger(__name__)

__all__ = [
    "DataSortingStage",
    "ExtraColumnRemovalStage",
    "FieldExtractionStage",
    "GeneBedCreationStage",
    "GenotypeReplacementStage",
    "MultiAllelicSplitStage",
    "PCAComputationStage",
    "ParallelCompleteProcessingStage",
    "PhenotypeIntegrationStage",
    "SnpSiftFilterStage",
    "TranscriptFilterStage",
    "VariantExtractionStage",
    "validate_annotation_consequence_filtering",
]


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
            regions_bed = context.config.get("regions_bed")

            if regions_bed and normalized_genes.lower().strip() == "all":
                context.gene_bed_file = validate_file_exists(regions_bed, self.name)
                logger.info(
                    "Using regions BED directly because gene selection is 'all': %s",
                    context.gene_bed_file,
                )
                return context

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
                # Check if it's a gene not found error before the broader tool-not-found case.
                if "not found in database" in str(e):
                    snpeff_reason = getattr(e, "stderr", None) or str(e)
                    raise FileFormatError(
                        normalized_genes,
                        (
                            f"valid gene names for {reference}; "
                            f"snpEff genes2bed reported: {snpeff_reason.strip()}"
                        ),
                        self.name,
                    ) from e
                # Check if it's a tool not found error
                if "snpEff" in str(e) and ("not found" in str(e) or "No such file" in str(e)):
                    raise ToolNotFoundError("snpEff", self.name) from e
                raise

            context.gene_bed_file = Path(bed_file)
            logger.info(f"Created BED file: {bed_file}")

            # Region restriction: intersect gene BED with restriction BED
            if regions_bed:
                context.gene_bed_file = self._intersect_with_restriction_bed(
                    context.gene_bed_file, regions_bed, context
                )
                logger.info(f"Gene BED restricted to regions: {context.gene_bed_file}")

            return context

    def _read_chromosomes(self, bed_file: Path) -> set[str]:
        """Read chromosome names from the first column of a BED file."""
        chromosomes: set[str] = set()
        with open(bed_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                chrom = line.split("\t")[0]
                if chrom:
                    chromosomes.add(chrom)
        return chromosomes

    def _count_regions(self, bed_file: Path) -> int:
        """Count non-comment, non-empty lines in a BED file."""
        count = 0
        with open(bed_file) as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith(("#", "track", "browser")):
                    continue
                count += 1
        return count

    def _intersect_with_restriction_bed(
        self, gene_bed: Path, restriction_bed: str, context: PipelineContext
    ) -> Path:
        """Intersect gene BED with a restriction BED file (e.g., capture kit)."""
        restriction_path = Path(restriction_bed)
        if not restriction_path.is_file():
            raise FileNotFoundError(f"Restriction BED file not found: {restriction_bed}")

        # Detect chr-prefix mismatch
        gene_chroms = self._read_chromosomes(gene_bed)
        restriction_chroms = self._read_chromosomes(restriction_path)

        if gene_chroms and restriction_chroms:
            gene_has_chr = any(c.startswith("chr") for c in gene_chroms)
            restriction_has_chr = any(c.startswith("chr") for c in restriction_chroms)
            if gene_has_chr != restriction_has_chr:
                raise ValueError(
                    "Chromosome naming mismatch between gene BED and restriction BED: "
                    f"gene BED {'has' if gene_has_chr else 'lacks'} 'chr' prefix "
                    f"(e.g., {next(iter(gene_chroms))}), "
                    f"restriction BED {'has' if restriction_has_chr else 'lacks'} 'chr' prefix "
                    f"(e.g., {next(iter(restriction_chroms))}). "
                    "Ensure both files use the same chromosome naming convention."
                )

        # Count original regions for logging
        orig_count = self._count_regions(gene_bed)

        # Run bedtools intersect
        out_path = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.restricted.bed"
        )
        with open(out_path, "w") as out_f:
            subprocess.run(
                ["bedtools", "intersect", "-a", str(gene_bed), "-b", str(restriction_bed)],
                stdout=out_f,
                check=True,
            )

        # Check result is non-empty
        if out_path.stat().st_size == 0:
            raise ValueError(
                f"Intersection of gene BED with restriction BED produced no regions. "
                f"Gene BED: {gene_bed}, Restriction BED: {restriction_bed}. "
                "Check that both files cover overlapping genomic regions."
            )

        new_count = self._count_regions(out_path)
        logger.info(
            f"Region restriction: {new_count} of {orig_count} gene regions retained "
            f"after intersection with {restriction_bed}"
        )

        return out_path

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


class MultiAllelicSplitStage(Stage):
    """Split multi-allelic SNPeff annotations into one per line."""

    def __init__(self, splitting_mode: str | None = None):
        """Initialize with the requested split timing mode."""
        super().__init__()
        self.splitting_mode = splitting_mode

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
        return {"variant_extraction"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        if self.splitting_mode == "after_filters":
            return {"snpsift_filtering", "transcript_filtering"}
        return set()

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return False

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Split SNPeff annotations if requested."""
        splitting_mode = self.splitting_mode or context.config.get("snpeff_splitting_mode")
        if not splitting_mode:
            logger.debug("SNPeff transcript splitting not requested")
            return context

        if splitting_mode == "after_filters":
            input_vcf = (
                getattr(context, "filtered_vcf", None)
                or getattr(context, "transcript_filtered_vcf", None)
                or getattr(context, "extracted_vcf", None)
                or context.data
            )
        else:
            input_vcf = getattr(context, "extracted_vcf", None) or context.data
        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.split_annotations.vcf.gz"
        )

        logger.info("Splitting SNPeff annotations by transcript")
        split_snpeff_annotations(str(input_vcf), str(output_vcf))

        context.split_annotations_vcf = output_vcf  # type: ignore[attr-defined]
        # Only update context.data if we haven't extracted to TSV yet
        if not getattr(context, "extracted_tsv", None):
            context.data = output_vcf
        return context


def validate_annotation_consequence_filtering(config: dict[str, Any]) -> None:
    """Validate that gene burden/association with consequence filters uses before_filters."""
    if config.get("snpeff_splitting_mode") == "before_filters":
        return

    analysis_requested = bool(
        config.get("perform_gene_burden") or config.get("perform_association")
    )
    if not analysis_requested:
        return

    filter_expr = config.get("filter") or config.get("filters") or ""
    has_consequence_filter = isinstance(filter_expr, str) and any(
        token in filter_expr for token in ("ANN[ANY]", "LOF[*]", "NMD[*]")
    )

    raw_presets = config.get("presets", [])
    presets = list(raw_presets) if isinstance(raw_presets, (list, tuple)) else [raw_presets]
    single_preset = config.get("preset")
    if single_preset:
        presets.extend(
            single_preset if isinstance(single_preset, (list, tuple)) else [single_preset]
        )

    consequence_presets = {"coding", "high", "missense", "high_or_lof_or_nmd", "damaging"}
    has_consequence_preset = any(p in consequence_presets for p in presets)

    if has_consequence_filter or has_consequence_preset:
        raise ValueError(
            "Gene burden and association analyses with consequence filtering "
            "(ANN[ANY], LOF[*], NMD[*]) require row-level qualification to prevent "
            "cross-gene misattribution. Please use '--split-snpeff-lines before_filters'."
        )


class SnpSiftFilterStage(Stage):
    """Apply SnpSift filter expressions."""

    def __init__(self, split_before_filtering: bool = False):
        """Initialize with whether split annotations should feed filtering."""
        super().__init__()
        self.split_before_filtering = split_before_filtering

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
        # Run after variant extraction and transcript filtering
        deps = {"variant_extraction", "transcript_filtering"}
        if self.split_before_filtering:
            deps.add("multiallelic_split")
        return deps

    def _split_before_filter(self, context: PipelineContext) -> bool:
        """Check if we should split before filtering."""
        # Check if snpeff_splitting_mode is set to "before_filters"
        splitting_mode = context.config.get("snpeff_splitting_mode", None)
        return bool(splitting_mode == "before_filters")

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Apply SnpSift filters."""
        validate_annotation_consequence_filtering(context.config)

        # Check both "filter" and "filters" keys (different parts use different names)
        filter_expr = context.config.get("filter") or context.config.get("filters")

        if not filter_expr or context.config.get("late_filtering"):
            logger.debug("No SnpSift filter or using late filtering")
            return context

        # Get input VCF from context (prefer transcript-filtered VCF when available)
        input_vcf = (
            getattr(context, "transcript_filtered_vcf", None)
            or (
                getattr(context, "split_annotations_vcf", None)
                if self.split_before_filtering
                else None
            )
            or getattr(context, "extracted_vcf", None)
            or context.data
        )
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
        filter_config = {
            "threads": context.config.get("threads", 1),
            "filter_exception_file": context.config.get("filter_exception_file"),
        }
        apply_snpsift_filter(str(input_vcf), filter_expr, filter_config, str(output_vcf))

        context.filtered_vcf = output_vcf
        # Only update context.data if we haven't extracted to TSV yet
        if not getattr(context, "extracted_tsv", None):
            context.data = output_vcf

        return context


class TranscriptFilterStage(Stage):
    """Filter SnpEff ANN annotations to requested transcript IDs."""

    def __init__(self, split_before_filtering: bool = False):
        """Initialize with whether split annotations should feed transcript filtering."""
        super().__init__()
        self.split_before_filtering = split_before_filtering

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "transcript_filtering"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Filter ANN annotations by transcript ID"

    @property
    def dependencies(self) -> set[str]:
        """Return the set of stage names this stage depends on."""
        return {"variant_extraction"}

    @property
    def soft_dependencies(self) -> set[str]:
        """Return the set of stage names that should run before if present."""
        if self.split_before_filtering:
            return {"multiallelic_split"}
        return set()

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Filter VCF ANN entries to requested transcript IDs."""
        transcript_ids = load_transcript_ids(
            context.config.get("transcript_list"),
            context.config.get("transcript_file"),
        )
        if not transcript_ids:
            logger.debug("No transcript filtering requested")
            return context

        input_vcf = (
            getattr(context, "split_annotations_vcf", None)
            or getattr(context, "extracted_vcf", None)
            or context.data
        )
        if not input_vcf:
            raise ValueError("No input VCF available for transcript filtering")

        output_vcf = context.workspace.get_intermediate_path(
            f"{context.workspace.base_name}.transcript_filtered.vcf.gz"
        )
        filter_vcf_to_transcripts(
            str(input_vcf),
            str(output_vcf),
            transcript_ids,
            index_output=False,
            filter_exception_file=context.config.get("filter_exception_file"),
        )

        context.transcript_filtered_vcf = output_vcf  # type: ignore[attr-defined]
        if not getattr(context, "extracted_tsv", None):
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
        return {"snpsift_filtering", "transcript_filtering", "multiallelic_split"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Extract fields to TSV format."""
        validate_annotation_consequence_filtering(context.config)
        # Use the most recent VCF in the processing chain
        if context.config.get("snpeff_splitting_mode") == "after_filters":
            input_vcf = (
                getattr(context, "split_annotations_vcf", None)
                or getattr(context, "filtered_vcf", None)
                or getattr(context, "transcript_filtered_vcf", None)
                or getattr(context, "extracted_vcf", None)
                or context.data
            )
        else:
            input_vcf = (
                getattr(context, "filtered_vcf", None)
                or getattr(context, "transcript_filtered_vcf", None)
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

        if context.config.get("filter_exception_file"):
            fields_str = ensure_fields_in_extract(fields_str, ["RAW_FILTER"])

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
    def parallel_safe(self) -> bool:
        return False

    @property
    def estimated_runtime(self) -> float:
        return 60.0

    def _validate_chunk_tsv(
        self,
        chunk_tsv: Path,
        fallback_on_error: bool = True,
        expected_header: str | None = None,
    ) -> bool:
        """Validate that a chunk TSV file is complete and valid."""
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
                    raw_line = f.readline()
                    first_line = (
                        raw_line.decode() if isinstance(raw_line, bytes) else str(raw_line)
                    ).strip()
                    # Check if it looks like a TSV header
                    if not first_line or "\t" not in first_line:
                        logger.debug(f"Chunk TSV {chunk_tsv} has invalid header format")
                        return False
                    if expected_header is not None and first_line != expected_header.strip():
                        return False
            except Exception as e:
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

        # Validate the merged TSV exists and is valid before mutating context.
        if not merged_tsv.exists() or not self._validate_chunk_tsv(merged_tsv):
            raise RuntimeError(
                "Cannot restore parallel_complete_processing from checkpoint: "
                f"merged TSV missing or invalid: {merged_tsv}"
            )
        logger.info(f"Restored valid merged TSV: {merged_tsv}")

        # Restore context state that would have been set
        context.extracted_tsv = merged_tsv
        context.data = merged_tsv
        context.config["parallel_vcf_processing_complete"] = True
        context.config.pop("chunked_processing_complete", None)

        # Mark all the stages this stage would have completed
        context.mark_complete("variant_extraction")
        context.mark_complete("snpsift_filtering")
        context.mark_complete("field_extraction")

        logger.info(
            f"Restored data path to {merged_tsv} and marked constituent stages as complete "
            f"after checkpoint skip"
        )
        return context

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Process variants in parallel with complete pipeline per chunk."""
        validate_annotation_consequence_filtering(context.config)
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
        context.mark_complete("snpsift_filtering")
        context.mark_complete("field_extraction")

        # Mark parallel VCF chunk processing as complete. DataFrame chunked
        # analysis has a separate completion flag.
        context.config["parallel_vcf_processing_complete"] = True
        context.config.pop("chunked_processing_complete", None)

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

        working_vcf = chunk_vcf
        split_mode = config.get("snpeff_splitting_mode")
        if split_mode and split_mode != "after_filters":
            split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
            split_snpeff_annotations(str(working_vcf), str(split_vcf))
            working_vcf = split_vcf

        transcript_ids = config.get("transcript_ids")
        if transcript_ids is None:
            transcript_ids = load_transcript_ids(
                config.get("transcript_list"),
                config.get("transcript_file"),
            )
        if transcript_ids:
            transcript_filtered_vcf = intermediate_dir / f"{chunk_base}.transcript_filtered.vcf.gz"
            filter_vcf_to_transcripts(
                str(working_vcf),
                str(transcript_filtered_vcf),
                transcript_ids,
                index_output=False,
                filter_exception_file=config.get("filter_exception_file"),
            )
            working_vcf = transcript_filtered_vcf

        # Step 2: Apply SnpSift filter (if not late filtering)
        filter_start = time.time()
        if config.get("late_filtering", False):
            filtered_vcf = working_vcf
            filter_time = 0.0
        else:
            filtered_vcf = intermediate_dir / f"{chunk_base}.filtered.vcf.gz"
            filter_expr = config.get("filter") or config.get("filters")
            if filter_expr:
                apply_snpsift_filter(
                    str(working_vcf),
                    filter_expr,
                    {
                        "threads": config.get("threads_per_chunk", 2),
                        "filter_exception_file": config.get("filter_exception_file"),
                    },
                    str(filtered_vcf),
                )
            else:
                filtered_vcf = working_vcf
        if split_mode == "after_filters":
            split_vcf = intermediate_dir / f"{chunk_base}.split_annotations.vcf.gz"
            split_snpeff_annotations(str(filtered_vcf), str(split_vcf))
            filtered_vcf = split_vcf
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

        if config.get("filter_exception_file"):
            fields_str = ensure_fields_in_extract(fields_str, ["RAW_FILTER"])

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
        # Handle empty bed_chunks case
        if not bed_chunks:
            logger.warning("No BED chunks to process - skipping parallel processing")
            return []

        vcf_file = context.config.get("vcf_file") or (str(context.data) if context.data else None)
        if not vcf_file:
            raise KeyError("vcf_file")
        threads = context.config.get("threads", 1)
        base_name = context.workspace.base_name
        intermediate_dir = context.workspace.intermediate_dir

        # Phase 1: Check for existing valid chunk TSVs
        chunks_to_process = []
        existing_chunk_tsvs = []

        use_compression = context.config.get("gzip_intermediates", True)
        tsv_suffix = ".extracted.tsv.gz" if use_compression else ".extracted.tsv"

        expected_header: str | None = None
        fields = context.config.get("fields_to_extract")
        if fields:
            cols = get_expected_extraction_columns(fields, getattr(context, "vcf_samples", None))
            expected_header = "\t".join(cols)

        for i, chunk_bed in enumerate(bed_chunks):
            chunk_base = f"{base_name}.chunk_{i}"
            expected_tsv = intermediate_dir / f"{chunk_base}{tsv_suffix}"

            if self._validate_chunk_tsv(expected_tsv, expected_header=expected_header):
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
                "snpeff_splitting_mode": context.config.get("snpeff_splitting_mode"),
                "filter": context.config.get("filter"),
                "filters": context.config.get("filters"),
                "filter_exception_file": context.config.get("filter_exception_file"),
                "fields_to_extract": context.config.get("fields_to_extract", ""),
                "append_extra_sample_fields": context.config.get(
                    "append_extra_sample_fields", False
                ),
                "extra_sample_fields": context.config.get("extra_sample_fields", []),
                "extract_fields_separator": context.config.get("extract_fields_separator", ","),
                "gzip_intermediates": use_compression,
                "vcf_samples": context.vcf_samples,
                "transcript_ids": load_transcript_ids(
                    context.config.get("transcript_list"),
                    context.config.get("transcript_file"),
                ),
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

        if not chunk_tsvs:
            fields = context.config.get("fields_to_extract", "")
            cols = get_expected_extraction_columns(fields, getattr(context, "vcf_samples", None))
            hdr = "\t".join(cols) + "\n"
            if use_compression:
                with gzip.open(output_tsv, "wt", compresslevel=1) as f:
                    f.write(hdr)
            else:
                output_tsv.write_text(hdr, encoding="utf-8")
            logger.info(f"Created schema-valid header-only TSV for empty chunks: {output_tsv}")
            return output_tsv

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
                    if first_header is None:
                        first_header = header.strip()
                        out_fh.write(header)
                    elif header.strip() != first_header:
                        raise ValueError(
                            f"Incompatible chunk header in {chunk_tsv}: "
                            f"expected '{first_header}', got '{header.strip()}'"
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

            # Remove chunks directory if empty
            if bed_chunks and bed_chunks[0].parent.name == "chunks":
                chunk_dir = bed_chunks[0].parent
                if chunk_dir.exists() and not any(chunk_dir.iterdir()):
                    chunk_dir.rmdir()

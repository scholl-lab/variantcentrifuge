"""
Processing stages for variant extraction and data transformation.

This module contains stages that handle the core data processing tasks:
- Gene BED file creation
- Variant extraction from VCF files
- Filtering with bcftools and SnpSift
- Field extraction to TSV format
- Phenotype integration
"""

import gzip  # noqa: F401
import logging
import os  # noqa: F401
import re
import shutil  # noqa: F401
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: F401
from pathlib import Path
from typing import Any

import pandas as pd  # noqa: F401

# smart_open.open exists in 7.x and 8.x; 8.x removed the old smart_open.smart_open alias.
from smart_open import open as smart_open  # noqa: F401

from ..extractor import extract_fields_bcftools, get_expected_extraction_columns  # noqa: F401
from ..filters import apply_snpsift_filter, extract_variants
from ..gene_bed import get_gene_bed, normalize_genes  # noqa: F401
from ..pipeline_core import PipelineContext, Stage
from ..pipeline_core.error_handling import (
    FileFormatError,
    ToolNotFoundError,
    graceful_error_handling,
    retry_on_failure,
    validate_file_exists,
)
from ..sample_subset import resolve_sample_subset
from ..transcript_filter import filter_vcf_to_transcripts, load_transcript_ids
from ..utils import ensure_fields_in_extract, split_bed_file  # noqa: F401
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


from .bed_stages import GeneBedCreationStage
from .parallel_processing_stage import ParallelCompleteProcessingStage


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
            validate_annotation_consequence_filtering(context.config)
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
            def extract() -> None:
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
    pattern = (
        r"\b(?:ANN|LOF|NMD)\s*\[\s*(?:ANY|\*|\d+)\s*\]"
        if config.get("snpeff_splitting_mode") == "after_filters"
        else r"\b(?:ANN\s*\[\s*(?:ANY|\*|[1-9]\d*)\s*\]|(?:LOF|NMD)\s*\[\s*(?:ANY|\*|\d+)\s*\])"
    )
    has_consequence_filter = isinstance(filter_expr, str) and bool(
        re.search(pattern, filter_expr, re.IGNORECASE)
    )

    selected_presets: list[str] = []
    single_preset = config.get("preset")
    if single_preset:
        if isinstance(single_preset, (list, tuple)):
            selected_presets.extend(p for p in single_preset if isinstance(p, str))
        elif isinstance(single_preset, str):
            selected_presets.append(single_preset)

    presets_list = config.get("presets_list")
    if presets_list:
        if isinstance(presets_list, (list, tuple)):
            selected_presets.extend(p for p in presets_list if isinstance(p, str))
        elif isinstance(presets_list, str):
            selected_presets.append(presets_list)

    raw_presets = config.get("presets")
    if raw_presets and not isinstance(raw_presets, dict):
        if isinstance(raw_presets, (list, tuple)):
            selected_presets.extend(p for p in raw_presets if isinstance(p, str))
        elif isinstance(raw_presets, str):
            selected_presets.append(raw_presets)

    consequence_presets = {"coding", "high", "missense", "high_or_lof_or_nmd", "damaging"}
    has_consequence_preset = any(p in consequence_presets for p in selected_presets)

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
            "max_memory_gb": context.config.get("max_memory_gb"),
            "snpsift_memory_gb": context.config.get("snpsift_memory_gb"),
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
        # Sample/phenotype stages fix the analysed sample set used by resolve_sample_subset.
        return {
            "snpsift_filtering",
            "transcript_filtering",
            "multiallelic_split",
            "sample_config_loading",
            "phenotype_case_control_assignment",
        }

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
            "bcftools_samples": resolve_sample_subset(context),
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

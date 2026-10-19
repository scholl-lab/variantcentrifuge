"""
Setup stages for configuration and data loading.

This module contains stages that load configuration files, phenotype data,
scoring configurations, and other setup tasks that can run in parallel.
"""

import logging
import os

from ..config import load_config, normalize_annotation_config
from ..phenotype import load_phenotypes
from ..pipeline_core import PipelineContext, Stage
from ..scoring import read_scoring_config
from .phenotype_case_control_stages import PhenotypeCaseControlAssignmentStage
from .sample_and_pedigree_stages import (
    PedigreeLoadingStage,
    SampleConfigLoadingStage,
    get_vcf_samples,
    read_pedigree,
)

logger = logging.getLogger(__name__)

ROW_LEVEL_CONSEQUENCE_TOKENS = ("ANN[ANY]", "LOF[*]", "NMD[*]")


def load_terms_from_file(file_path: str | None, logger: logging.Logger) -> list[str]:
    """
    Load terms (HPO terms, sample IDs, etc.) from a file, one per line.

    Parameters
    ----------
    file_path : str or None
        Path to a file containing one term per line.
    logger : logging.Logger
        Logger instance for error logging.

    Returns
    -------
    list of str
        A list of terms loaded from the file.

    Raises
    ------
    SystemExit
        If the file is missing or empty and a file_path was specified.
    """
    terms: list[str] = []
    if file_path:
        if not os.path.exists(file_path):
            logger.error(f"Required file not found: {file_path}")
            raise FileNotFoundError(f"Required file not found: {file_path}")
        with open(file_path, encoding="utf-8") as f:
            found_any = False
            for line in f:
                t = line.strip()
                if t:
                    found_any = True
                    terms.append(t)
            if not found_any:
                logger.error(f"File {file_path} is empty or invalid.")
                raise ValueError(f"File {file_path} is empty or invalid.")
    return terms


def _warn_if_record_level_consequence_mode(config: dict) -> None:
    """Warn when analysis options request row-level semantics from a record-level mode."""
    if config.get("snpeff_splitting_mode") == "before_filters":
        return

    transcript_filtering_requested = bool(
        config.get("transcript_list")
        or config.get("transcript_file")
        or config.get("transcript_ids")
    )
    analysis_requested = bool(
        config.get("perform_gene_burden") or config.get("perform_association")
    )
    filter_expr = config.get("filter") or config.get("filters") or ""
    consequence_filter_requested = isinstance(filter_expr, str) and any(
        token in filter_expr for token in ROW_LEVEL_CONSEQUENCE_TOKENS
    )

    if transcript_filtering_requested and analysis_requested and consequence_filter_requested:
        logger.warning(
            "Transcript-filtered gene burden/association with ANN[ANY], LOF[*], or "
            "NMD[*] consequence predicates should use --split-snpeff-lines before_filters "
            "for row-level consequence qualification. Without before_filters, SnpSift "
            "evaluates those predicates on the original multi-annotation record before "
            "transcript selection."
        )


class ConfigurationLoadingStage(Stage):
    """Load and merge configuration from file and CLI arguments."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "configuration_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Configuration is already present in the new CLI context."""
        return "skip"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load configuration and merge with CLI arguments"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load configuration and merge with CLI arguments."""
        logger.info("Loading configuration...")

        # Check if config was already passed as a dict (from main pipeline)
        if hasattr(context.args, "config") and isinstance(context.args.config, dict):
            # Config already fully merged by cli.py, just use it
            config = context.args.config
            logger.info("Using pre-merged configuration from main pipeline")
        else:
            # Load base configuration
            if hasattr(context.args, "config") and context.args.config:
                config = load_config(context.args.config)
                logger.info(f"Loaded configuration from {context.args.config}")
            else:
                config = load_config()  # Load default config
                logger.info("Loaded default configuration")

            # Merge CLI arguments (CLI takes precedence)
            cli_args = vars(context.args)
            for key, value in cli_args.items():
                if value is not None and key not in ["config", "fields"]:
                    # Skip internal arguments and fields (handled as fields_to_extract)
                    if key.startswith("_"):
                        continue
                    config[key] = value

            # Ensure phenotype arguments are properly set even if they come as None
            # This fixes the issue where --phenotype-file args don't make it to config
            phenotype_args = ["phenotype_file", "phenotype_sample_column", "phenotype_value_column"]
            for arg in phenotype_args:
                if hasattr(context.args, arg):
                    config[arg] = getattr(context.args, arg)
        # Merge any existing config from context (e.g., from tests)
        if context.config:
            for key, value in context.config.items():
                if key not in config:
                    config[key] = value
        # Set pipeline version - use actual package version
        from ..version import __version__

        config["pipeline_version"] = config.get("pipeline_version", __version__)
        # Convert fields_to_extract to extract for downstream stages
        if "extract" not in config and "fields_to_extract" in config:
            config["extract"] = config["fields_to_extract"].split()
            logger.debug(
                f"Converted fields_to_extract to extract list with {len(config['extract'])} fields"
            )

        # Parse phenotype arguments from CLI
        # This fixes the bug where phenotype-based case/control assignment doesn't work
        case_hpo_terms = []
        control_hpo_terms = []

        # Parse case phenotypes
        if config.get("case_phenotypes"):
            if isinstance(config["case_phenotypes"], str):
                # Parse comma-separated string
                case_hpo_terms = [
                    t.strip() for t in config["case_phenotypes"].split(",") if t.strip()
                ]
                logger.debug(f"Parsed {len(case_hpo_terms)} case phenotype terms from CLI")
            else:
                # Already a list
                case_hpo_terms = config["case_phenotypes"]
        # Load additional case phenotypes from file
        if config.get("case_phenotypes_file"):
            file_terms = load_terms_from_file(config["case_phenotypes_file"], logger)
            case_hpo_terms.extend(file_terms)
            logger.debug(f"Loaded {len(file_terms)} additional case phenotype terms from file")

        # Parse control phenotypes
        if config.get("control_phenotypes"):
            if isinstance(config["control_phenotypes"], str):
                # Parse comma-separated string
                control_hpo_terms = [
                    t.strip() for t in config["control_phenotypes"].split(",") if t.strip()
                ]
                logger.debug(f"Parsed {len(control_hpo_terms)} control phenotype terms from CLI")
            else:
                # Already a list
                control_hpo_terms = config["control_phenotypes"]

        # Load additional control phenotypes from file
        if config.get("control_phenotypes_file"):
            file_terms = load_terms_from_file(config["control_phenotypes_file"], logger)
            control_hpo_terms.extend(file_terms)
            logger.debug(f"Loaded {len(file_terms)} additional control phenotype terms from file")

        config["case_phenotypes"] = case_hpo_terms
        config["control_phenotypes"] = control_hpo_terms

        case_samples = []
        control_samples = []

        if config.get("case_samples"):
            if isinstance(config["case_samples"], str):
                case_samples = [s.strip() for s in config["case_samples"].split(",") if s.strip()]
                logger.debug(f"Parsed {len(case_samples)} case samples from CLI")
            else:
                case_samples = config["case_samples"]

        if config.get("control_samples"):
            if isinstance(config["control_samples"], str):
                control_samples = [
                    s.strip() for s in config["control_samples"].split(",") if s.strip()
                ]
                logger.debug(f"Parsed {len(control_samples)} control samples from CLI")
            else:
                control_samples = config["control_samples"]

        config["case_samples"] = case_samples
        config["control_samples"] = control_samples

        _warn_if_record_level_consequence_mode(config)
        from .processing_stages import validate_annotation_consequence_filtering

        validate_annotation_consequence_filtering(config)

        context.config = config

        logger.debug(f"Configuration loaded with {len(config)} settings")
        return context


class PhenotypeLoadingStage(Stage):
    """Load phenotype data from file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "phenotype_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Phenotype data can be restored by the checkpoint skip handler."""
        return "restore"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load phenotype data"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load phenotype data if provided."""
        phenotype_file = context.config.get("phenotype_file")
        phenotype_sample_column = context.config.get("phenotype_sample_column")
        phenotype_value_column = context.config.get("phenotype_value_column")

        logger.debug(f"Phenotype file: {phenotype_file}")
        logger.debug(f"Phenotype sample column: {phenotype_sample_column}")
        logger.debug(f"Phenotype value column: {phenotype_value_column}")
        logger.debug(f"Config keys: {list(context.config.keys())}")

        if phenotype_file and phenotype_sample_column and phenotype_value_column:
            logger.info(f"Loading phenotypes from {phenotype_file}")

            phenotypes = load_phenotypes(
                phenotype_file, phenotype_sample_column, phenotype_value_column
            )

            if not phenotypes:
                raise ValueError(
                    f"No phenotype data loaded from {phenotype_file}. "
                    "Check file formatting and column names."
                )

            context.phenotype_data = phenotypes
            context.config["phenotypes"] = phenotypes  # For compatibility
            logger.info(f"Loaded phenotypes for {len(phenotypes)} samples")
        else:
            logger.debug("No phenotype file specified, skipping phenotype loading")

        return context

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        Since this is the phenotype loading stage, we need to restore phenotype data
        that would have been loaded during normal execution.
        """
        phenotype_file = context.config.get("phenotype_file")
        phenotype_sample_column = context.config.get("phenotype_sample_column")
        phenotype_value_column = context.config.get("phenotype_value_column")

        if phenotype_file and phenotype_sample_column and phenotype_value_column:
            logger.debug(f"Restoring phenotypes from {phenotype_file} (checkpoint skip)")

            phenotypes = load_phenotypes(
                phenotype_file, phenotype_sample_column, phenotype_value_column
            )

            if phenotypes:
                context.phenotype_data = phenotypes
                context.config["phenotypes"] = phenotypes  # For compatibility
                logger.debug(f"Restored phenotypes for {len(phenotypes)} samples (checkpoint skip)")
            else:
                raise ValueError(
                    f"No phenotype data loaded from {phenotype_file}. "
                    "Check file formatting and column names."
                )
        else:
            logger.debug("No phenotype configuration found during checkpoint skip")

        return context


class ScoringConfigLoadingStage(Stage):
    """Load variant scoring configuration."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "scoring_config_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Scoring configuration can be restored by the checkpoint skip handler."""
        return "restore"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load scoring configuration"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load scoring configuration if provided."""
        scoring_config_path = context.config.get("scoring_config_path")

        if scoring_config_path:
            logger.info(f"Loading scoring configuration from {scoring_config_path}")

            try:
                scoring_config = read_scoring_config(scoring_config_path)
                context.scoring_config = scoring_config
                logger.info(
                    f"Successfully loaded scoring configuration with "
                    f"{len(scoring_config.get('formulas', {}))} formulas"
                )
            except Exception as e:
                raise ValueError(f"Failed to load scoring configuration: {e}") from e
        else:
            logger.debug("No scoring configuration specified")

        return context

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        Since this is the scoring config loading stage, we need to restore scoring
        configuration that would have been loaded during normal execution.
        """
        scoring_config_path = context.config.get("scoring_config_path")

        if scoring_config_path:
            logger.debug(
                f"Restoring scoring configuration from {scoring_config_path} (checkpoint skip)"
            )

            try:
                scoring_config = read_scoring_config(scoring_config_path)
                context.scoring_config = scoring_config
                logger.debug(
                    f"Restored scoring configuration with "
                    f"{len(scoring_config.get('formulas', {}))} formulas (checkpoint skip)"
                )
            except Exception as e:
                raise ValueError(f"Failed to restore scoring configuration: {e}") from e
        else:
            logger.debug("No scoring configuration specified during checkpoint skip")

        return context


class AnnotationConfigLoadingStage(Stage):
    """Load custom annotation configurations (BED files, gene lists, JSON)."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "annotation_config_loading"

    @property
    def checkpoint_resume_policy(self) -> str:
        """Annotation configuration can be restored by the checkpoint skip handler."""
        return "restore"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load annotation configurations"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _sync_args_to_config(self, context: PipelineContext) -> None:
        """Sync annotation CLI arguments into context.config if not already present."""
        if hasattr(context, "args") and context.args:
            keys = [
                "annotate_bed",
                "annotate_gene_list",
                "annotate_json",
                "annotate_json_genes",
            ]
            for key in keys:
                val = getattr(context.args, key, None)
                if val:
                    context.config.setdefault(key, val)
        normalize_annotation_config(context.config)

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load all annotation configurations."""
        annotations = {}
        self._sync_args_to_config(context)

        # Load BED file annotations
        bed_files = context.config.get("annotate_bed_files", [])
        if bed_files:
            if isinstance(bed_files, str):
                bed_files = [bed_files]
            annotations["bed_files"] = bed_files
            logger.info(f"Configured {len(bed_files)} BED file annotations")

        # Load gene list annotations
        gene_lists = context.config.get("annotate_gene_lists", [])
        if gene_lists:
            if isinstance(gene_lists, str):
                gene_lists = [gene_lists]
            annotations["gene_lists"] = gene_lists
            logger.info(f"Configured {len(gene_lists)} gene list annotations")

        # Load JSON gene annotations
        json_genes = context.config.get("annotate_json_genes")
        json_mapping = context.config.get("json_gene_mapping")
        if json_genes:
            if isinstance(json_genes, str):
                json_genes = [json_genes]
            annotations["json_genes"] = json_genes
            if json_mapping:
                annotations["json_mapping"] = json_mapping
            logger.info(f"Configured JSON gene annotations from {json_genes}")

        if annotations:
            context.annotation_configs = annotations

        return context

    def _handle_checkpoint_skip(self, context: PipelineContext) -> PipelineContext:
        """Handle the case where this stage is skipped by checkpoint system.

        Since this is the annotation config loading stage, we need to restore annotation
        configurations that would have been loaded during normal execution.
        """
        annotations = {}
        self._sync_args_to_config(context)

        # Restore BED file annotations
        bed_files = context.config.get("annotate_bed_files", [])
        if bed_files:
            if isinstance(bed_files, str):
                bed_files = [bed_files]
            annotations["bed_files"] = bed_files
            logger.debug(f"Restored {len(bed_files)} BED file annotations (checkpoint skip)")

        # Restore gene list annotations
        gene_lists = context.config.get("annotate_gene_lists", [])
        if gene_lists:
            if isinstance(gene_lists, str):
                gene_lists = [gene_lists]
            annotations["gene_lists"] = gene_lists
            logger.debug(f"Restored {len(gene_lists)} gene list annotations (checkpoint skip)")

        # Restore JSON gene annotations
        json_genes = context.config.get("annotate_json_genes")
        json_mapping = context.config.get("json_gene_mapping")
        if json_genes:
            if isinstance(json_genes, str):
                json_genes = [json_genes]
            annotations["json_genes"] = json_genes
            if json_mapping:
                annotations["json_mapping"] = json_mapping
            logger.debug(f"Restored JSON gene annotations from {json_genes} (checkpoint skip)")

        if annotations:
            context.annotation_configs = annotations
            logger.debug(
                f"Restored annotation configurations with {len(annotations)} types "
                f"(checkpoint skip)"
            )

        return context


__all__ = [
    "AnnotationConfigLoadingStage",
    "ConfigurationLoadingStage",
    "PedigreeLoadingStage",
    "PhenotypeCaseControlAssignmentStage",
    "PhenotypeLoadingStage",
    "SampleConfigLoadingStage",
    "ScoringConfigLoadingStage",
    "get_vcf_samples",
    "load_terms_from_file",
    "read_pedigree",
]

"""
Setup stages for configuration and data loading.

This module contains stages that load configuration files, phenotype data,
scoring configurations, and other setup tasks that can run in parallel.
"""

import logging

from ..config import load_config
from ..helpers import get_vcf_samples
from ..ped_reader import read_pedigree
from ..phenotype import load_phenotypes
from ..pipeline_core import PipelineContext, Stage
from ..scoring import read_scoring_config

logger = logging.getLogger(__name__)


class ConfigurationLoadingStage(Stage):
    """Load and merge configuration from file and CLI arguments."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "configuration_loading"

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

        # Merge any existing config from context (e.g., from tests)
        if context.config:
            for key, value in context.config.items():
                if key not in config:
                    config[key] = value

        # Set pipeline version
        config["pipeline_version"] = config.get("pipeline_version", "1.0.0")

        # Convert fields_to_extract to extract for downstream stages
        if "extract" not in config and "fields_to_extract" in config:
            config["extract"] = config["fields_to_extract"].split()
            logger.debug(
                f"Converted fields_to_extract to extract list with {len(config['extract'])} fields"
            )

        # Update context
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


class ScoringConfigLoadingStage(Stage):
    """Load variant scoring configuration."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "scoring_config_loading"

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
                raise ValueError(f"Failed to load scoring configuration: {e}")
        else:
            logger.debug("No scoring configuration specified")

        return context


class PedigreeLoadingStage(Stage):
    """Load pedigree data from PED file."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "pedigree_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load pedigree data"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load pedigree data if provided."""
        ped_file = context.config.get("ped_file") or context.config.get("ped")

        if ped_file:
            logger.info(f"Loading pedigree from {ped_file}")

            try:
                pedigree_data = read_pedigree(ped_file)
                context.pedigree_data = pedigree_data

                # Log pedigree summary
                n_samples = len(pedigree_data)
                n_affected = sum(
                    1 for s in pedigree_data.values() if s.get("affected_status") == "2"
                )
                logger.info(f"Loaded pedigree with {n_samples} samples ({n_affected} affected)")
            except Exception as e:
                raise ValueError(f"Failed to load pedigree file: {e}")
        else:
            logger.debug("No pedigree file specified")

        return context


class AnnotationConfigLoadingStage(Stage):
    """Load custom annotation configurations (BED files, gene lists, JSON)."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "annotation_config_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load annotation configurations"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load all annotation configurations."""
        annotations = {}

        # Load BED file annotations
        bed_files = context.config.get("annotate_bed", [])
        if bed_files:
            if isinstance(bed_files, str):
                bed_files = [bed_files]
            annotations["bed_files"] = bed_files
            logger.info(f"Configured {len(bed_files)} BED file annotations")

        # Load gene list annotations
        gene_lists = context.config.get("annotate_gene_list", [])
        if gene_lists:
            if isinstance(gene_lists, str):
                gene_lists = [gene_lists]
            annotations["gene_lists"] = gene_lists
            logger.info(f"Configured {len(gene_lists)} gene list annotations")

        # Load JSON gene annotations
        json_genes = context.config.get("annotate_json_genes")
        json_mapping = context.config.get("json_gene_mapping")
        if json_genes and json_mapping:
            annotations["json_genes"] = json_genes
            annotations["json_mapping"] = json_mapping
            logger.info(f"Configured JSON gene annotations from {json_genes}")

        if annotations:
            context.annotation_configs = annotations

        return context


class SampleConfigLoadingStage(Stage):
    """Load sample configurations and VCF sample names."""

    @property
    def name(self) -> str:
        """Return the stage name."""
        return "sample_config_loading"

    @property
    def description(self) -> str:
        """Return a description of what this stage does."""
        return "Load sample configurations"

    @property
    def parallel_safe(self) -> bool:
        """Return whether this stage can run in parallel with others."""
        return True

    @property
    def dependencies(self) -> set:
        """Return the set of stage names this stage depends on."""
        # Needs config to know VCF file path
        return {"configuration_loading"}

    def _process(self, context: PipelineContext) -> PipelineContext:
        """Load VCF samples and any sample-specific configurations."""
        vcf_file = context.config.get("vcf_file")

        if not vcf_file:
            raise ValueError("No VCF file specified in configuration")

        # Get VCF samples
        logger.info(f"Reading sample names from {vcf_file}")
        vcf_samples = get_vcf_samples(vcf_file)

        # Apply sample substring removal if configured
        remove_substring = context.config.get("remove_sample_substring")
        if remove_substring and remove_substring.strip():
            logger.info(f"Removing substring '{remove_substring}' from sample names")
            original_samples = vcf_samples[:]  # Create a copy for logging
            vcf_samples = [s.replace(remove_substring, "") for s in vcf_samples]

            # Log changes for debugging
            changes = []
            for orig, new in zip(original_samples, vcf_samples):
                if orig != new:
                    changes.append(f"{orig} -> {new}")

            if changes:
                logger.info(f"Applied substring removal to {len(changes)} samples")
                logger.debug(f"Sample name changes: {changes[:5]}...")  # Log first 5 changes
            else:
                logger.info(f"No samples contained substring '{remove_substring}'")

        context.vcf_samples = vcf_samples
        logger.info(f"Found {len(vcf_samples)} samples in VCF")

        # Load case/control sample lists if provided
        case_samples_file = context.config.get("case_samples_file")
        control_samples_file = context.config.get("control_samples_file")

        if case_samples_file:
            with open(case_samples_file, "r") as f:
                case_samples = [line.strip() for line in f if line.strip()]
            context.config["case_samples"] = case_samples
            logger.info(f"Loaded {len(case_samples)} case samples")
            logger.debug(f"Case samples: {case_samples[:10]}...")  # Log first 10 for debugging

        if control_samples_file:
            with open(control_samples_file, "r") as f:
                control_samples = [line.strip() for line in f if line.strip()]
            context.config["control_samples"] = control_samples
            logger.info(f"Loaded {len(control_samples)} control samples")
            logger.debug(
                f"Control samples: {control_samples[:10]}..."
            )  # Log first 10 for debugging

        return context

"""
Refactored pipeline using the new Stage architecture.

This module demonstrates how to use the refactored pipeline stages
to replace the monolithic pipeline.py functionality.
"""

import argparse
import logging
from pathlib import Path
from typing import List

from .config import load_config
from .gene_bed import normalize_genes
from .pipeline_core.context import PipelineContext
from .pipeline_core.runner import PipelineRunner
from .pipeline_core.workspace import Workspace
from .stages.analysis_stages import (
    ChunkedAnalysisStage,
    CustomAnnotationStage,
    DataFrameLoadingStage,
    GeneBurdenAnalysisStage,
    InheritanceAnalysisStage,
    StatisticsGenerationStage,
    VariantAnalysisStage,
    VariantScoringStage,
)
from .stages.output_stages import (
    ArchiveCreationStage,
    ExcelReportStage,
    FinalFilteringStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    PseudonymizationStage,
    TSVOutputStage,
    VariantIdentifierStage,
)
from .stages.processing_stages import (
    DataSortingStage,
    ExtraColumnRemovalStage,
    FieldExtractionStage,
    GeneBedCreationStage,
    GenotypeReplacementStage,
    MultiAllelicSplitStage,
    ParallelCompleteProcessingStage,
    PhenotypeIntegrationStage,
    SnpSiftFilterStage,
    VariantExtractionStage,
)
from .stages.setup_stages import (
    AnnotationConfigLoadingStage,
    ConfigurationLoadingStage,
    PedigreeLoadingStage,
    PhenotypeCaseControlAssignmentStage,
    PhenotypeLoadingStage,
    SampleConfigLoadingStage,
    ScoringConfigLoadingStage,
)
from .utils import compute_base_name

logger = logging.getLogger(__name__)


def check_scoring_requires_inheritance(args: argparse.Namespace, config: dict) -> bool:
    """Check if scoring configuration requires inheritance analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
    config : dict
        Configuration dictionary

    Returns
    -------
    bool
        True if inheritance analysis is required for scoring
    """
    # Check if there's a scoring config
    scoring_config_path = config.get("scoring_config_path") or getattr(
        args, "scoring_config_path", None
    )
    if not scoring_config_path:
        return False

    try:
        # Import here to avoid circular imports
        from .scoring import read_scoring_config

        # Load scoring config to check for inheritance dependencies
        scoring_config = read_scoring_config(scoring_config_path)

        # Check if any scoring formulas reference inheritance variables
        if scoring_config and "formulas" in scoring_config:
            for formula_dict in scoring_config["formulas"]:
                for formula_name, formula_expr in formula_dict.items():
                    if isinstance(formula_expr, str):
                        # Check if formula references common inheritance variables
                        if any(
                            var in formula_expr
                            for var in [
                                "pattern",
                                "details",
                                "Inheritance_Pattern",
                                "Inheritance_Details",
                            ]
                        ):
                            logger.info(
                                f"Scoring formula '{formula_name}' requires inheritance analysis"
                            )
                            return True

        # Check variable assignments for inheritance column dependencies
        if scoring_config and "variables" in scoring_config:
            variables = scoring_config["variables"]
            if any(col in variables for col in ["Inheritance_Pattern", "Inheritance_Details"]):
                logger.info("Scoring configuration requires inheritance columns")
                return True

    except Exception as e:
        logger.warning(f"Could not check scoring config for inheritance dependencies: {e}")

    return False


def build_pipeline_stages(args: argparse.Namespace) -> List:
    """Build the list of stages based on configuration.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments

    Returns
    -------
    List[Stage]
        List of stages to execute
    """
    stages = []

    # Always load configuration
    stages.append(ConfigurationLoadingStage())

    # Get config to check for scoring_config_path
    # If args.config is a dict, use it; otherwise check individual args
    config = args.config if hasattr(args, "config") and isinstance(args.config, dict) else {}

    # Setup stages (can run in parallel)
    # Always add phenotype loading stage - it will skip if no phenotype file is provided
    stages.append(PhenotypeLoadingStage())

    # Check for scoring config in merged config or args
    if config.get("scoring_config_path") or (
        hasattr(args, "scoring_config_path") and args.scoring_config_path
    ):
        stages.append(ScoringConfigLoadingStage())

    if hasattr(args, "ped_file") and args.ped_file:
        stages.append(PedigreeLoadingStage())

    if any(
        [
            hasattr(args, "annotate_bed") and args.annotate_bed,
            hasattr(args, "annotate_gene_list") and args.annotate_gene_list,
            hasattr(args, "annotate_json_genes") and args.annotate_json_genes,
        ]
    ):
        stages.append(AnnotationConfigLoadingStage())

    # Always load samples
    stages.append(SampleConfigLoadingStage())

    # Add phenotype-based case/control assignment stage
    # This must run after both phenotype loading and sample config loading
    stages.append(PhenotypeCaseControlAssignmentStage())

    # Processing stages
    stages.append(GeneBedCreationStage())

    # Choose parallel or single-threaded processing
    threads = getattr(args, "threads", 1)
    if threads > 1:
        # Use the new complete parallel processing stage that runs the entire
        # pipeline (extraction -> filtering -> field extraction) for each chunk
        stages.append(ParallelCompleteProcessingStage())
    else:
        # Single-threaded processing uses individual stages
        stages.append(VariantExtractionStage())

        # Optional processing stages
        if hasattr(args, "snpeff_split_by_transcript") and args.snpeff_split_by_transcript:
            stages.append(MultiAllelicSplitStage())

        # Add SnpSift filtering stage (it will check config for filters at runtime)
        if not getattr(args, "late_filtering", False):
            stages.append(SnpSiftFilterStage())

        # Field extraction
        stages.append(FieldExtractionStage())

        # Data sorting (for efficient chunked processing)
        stages.append(DataSortingStage())

    # Optional transformations
    # Default behavior: do genotype replacement unless --no-replacement is specified
    if not getattr(args, "no_replacement", False):
        stages.append(GenotypeReplacementStage())

    # Always add phenotype integration stage - it will skip if no phenotype data is loaded
    stages.append(PhenotypeIntegrationStage())

    if hasattr(args, "extra_columns_to_remove") and args.extra_columns_to_remove:
        stages.append(ExtraColumnRemovalStage())

    # Analysis stages
    stages.append(DataFrameLoadingStage())
    # Add chunked analysis stage that will handle all analysis steps if chunked processing is needed
    # This stage will automatically activate when use_chunked_processing is set to True
    stages.append(ChunkedAnalysisStage())

    if any(
        [
            hasattr(args, "annotate_bed") and args.annotate_bed,
            hasattr(args, "annotate_gene_list") and args.annotate_gene_list,
            hasattr(args, "annotate_json_genes") and args.annotate_json_genes,
        ]
    ):
        stages.append(CustomAnnotationStage())

    # Check both args and config for inheritance settings
    # Debug logging to understand the issue
    has_ped_file_arg = hasattr(args, "ped_file") and args.ped_file
    has_inheritance_mode_arg = hasattr(args, "inheritance_mode") and args.inheritance_mode
    has_calculate_inheritance_config = config.get("calculate_inheritance", False)
    has_ped_file_config = config.get("ped_file")
    has_inheritance_mode_config = config.get("inheritance_mode")
    requires_inheritance_for_scoring = check_scoring_requires_inheritance(args, config)

    logger.debug("Inheritance analysis conditions:")
    logger.debug(f"  has_ped_file_arg: {has_ped_file_arg}")
    logger.debug(f"  has_inheritance_mode_arg: {has_inheritance_mode_arg}")
    logger.debug(f"  has_calculate_inheritance_config: {has_calculate_inheritance_config}")
    logger.debug(f"  has_ped_file_config: {has_ped_file_config}")
    logger.debug(f"  has_inheritance_mode_config: {has_inheritance_mode_config}")
    logger.debug(f"  requires_inheritance_for_scoring: {requires_inheritance_for_scoring}")

    should_calculate_inheritance = (
        has_ped_file_arg
        or has_inheritance_mode_arg
        or has_calculate_inheritance_config
        or has_ped_file_config
        or has_inheritance_mode_config
        or requires_inheritance_for_scoring
    )

    logger.info(f"Should calculate inheritance: {should_calculate_inheritance}")

    if should_calculate_inheritance:
        stages.append(InheritanceAnalysisStage())

    # Always run variant analysis to add standard columns
    stages.append(VariantAnalysisStage())

    # Add VAR_ID after variant analysis since it creates a new DataFrame
    stages.append(VariantIdentifierStage())

    # Check for scoring config in merged config or args
    # IMPORTANT: Scoring must happen AFTER variant analysis and variant identifier
    # because it may depend on columns created by those stages
    if config.get("scoring_config_path") or (
        hasattr(args, "scoring_config_path") and args.scoring_config_path
    ):
        stages.append(VariantScoringStage())

    if not getattr(args, "no_stats", False):
        stages.append(StatisticsGenerationStage())

    if hasattr(args, "perform_gene_burden") and args.perform_gene_burden:
        stages.append(GeneBurdenAnalysisStage())

    # Check for final filtering - need to check both args attributes and config
    late_filtering = getattr(args, "late_filtering", False) or config.get("late_filtering", False)
    final_filter = getattr(args, "final_filter", None) or config.get("final_filter", None)

    logger.debug(
        f"Checking for final filtering - late_filtering: {late_filtering}, "
        f"final_filter: {final_filter}"
    )

    if late_filtering or final_filter:
        logger.debug("Adding FinalFilteringStage to pipeline")
        stages.append(FinalFilteringStage())

    if hasattr(args, "pseudonymize") and args.pseudonymize:
        stages.append(PseudonymizationStage())

    stages.append(TSVOutputStage())

    # Report generation - add individual stages to main pipeline
    # The PipelineRunner will handle dependencies and execution order correctly
    if getattr(args, "xlsx", False) or getattr(args, "excel", False):
        stages.append(ExcelReportStage())

    if getattr(args, "html_report", False):
        stages.append(HTMLReportStage())

    if getattr(args, "igv", False):
        stages.append(IGVReportStage())

    if not getattr(args, "no_metadata", False):
        stages.append(MetadataGenerationStage())

    if hasattr(args, "archive_results") and args.archive_results:
        stages.append(ArchiveCreationStage())

    return stages


def create_stages_from_config(config: dict) -> List:
    """Create pipeline stages from a configuration dictionary.

    This function is used by CLI handlers that need to create stages
    without full argparse.Namespace objects.

    Parameters
    ----------
    config : dict
        Configuration dictionary with pipeline settings

    Returns
    -------
    List[Stage]
        List of stages for the current configuration
    """
    # Convert config dict to a minimal args namespace for compatibility
    args = argparse.Namespace()

    # Set essential attributes from config
    args.config = config
    args.scoring_config_path = config.get("scoring_config_path")
    args.ped_file = config.get("ped_file")
    args.annotate_bed = config.get("annotate_bed", [])
    args.annotate_gene_list = config.get("annotate_gene_list", [])
    args.annotate_json_genes = config.get("annotate_json_genes", [])
    args.chunks = config.get("chunks")
    args.late_filtering = config.get("late_filtering", False)
    args.genotype_filter = config.get("genotype_filter")
    args.final_filter = config.get("final_filter")
    args.pseudonymize = config.get("pseudonymize", False)
    args.xlsx = config.get("xlsx", False)
    args.html_report = config.get("html_report", False)
    args.igv_report = config.get("igv_report", False)
    args.archive_results = config.get("archive_results", False)

    # Use the existing build_pipeline_stages function
    return build_pipeline_stages(args)


def run_refactored_pipeline(args: argparse.Namespace) -> None:
    """Run the refactored pipeline with the new architecture.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
    """
    # Initialize stage registry for selective resume functionality
    from .stages.stage_registry import initialize_registry

    initialize_registry()

    # Setup logging
    log_level = getattr(args, "log_level", "INFO")
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Normalize genes first to compute base name
    gene_name = normalize_genes(
        getattr(args, "gene_name", "") or "", getattr(args, "gene_file", None), logger
    )

    # Compute base name like the original pipeline
    vcf_file = getattr(args, "vcf_file", "")
    base_name = compute_base_name(vcf_file, gene_name)

    # Create workspace
    output_dir = getattr(args, "output_dir", "output")
    workspace = Workspace(Path(output_dir), base_name)

    # Load initial config - it might be passed as args.config or already merged
    logger.debug(
        f"Checking args.config: hasattr={hasattr(args, 'config')}, "
        f"type={type(getattr(args, 'config', None)) if hasattr(args, 'config') else 'N/A'}"
    )
    if hasattr(args, "config") and isinstance(args.config, dict):
        # Config already loaded and passed as dict
        initial_config = args.config
        logger.debug(f"Using args.config as dict with {len(initial_config)} keys")
    elif hasattr(args, "config") and args.config:
        # Config file path
        initial_config = load_config(args.config)
        logger.debug(f"Loaded config from file: {args.config}")
    else:
        initial_config = {}
        logger.debug("Using empty initial config")

    # Ensure config has extract field from fields_to_extract
    if "extract" not in initial_config and "fields_to_extract" in initial_config:
        initial_config["extract"] = initial_config["fields_to_extract"].split()

    # Debug log to understand config state
    logger.debug(f"Initial config keys: {list(initial_config.keys())}")
    logger.debug(f"fields_to_extract in config: {'fields_to_extract' in initial_config}")
    logger.debug(f"extract in config: {'extract' in initial_config}")
    logger.debug(f"final_filter in config: {'final_filter' in initial_config}")
    if "fields_to_extract" in initial_config:
        logger.debug(f"fields_to_extract value: {initial_config['fields_to_extract'][:100]}...")
    if "extract" in initial_config:
        logger.debug(f"extract fields count: {len(initial_config['extract'])}")
    if "final_filter" in initial_config:
        logger.debug(f"final_filter value: {initial_config['final_filter']}")

    # Add normalized genes to config
    initial_config["normalized_genes"] = gene_name

    # Add gene_name and gene_file from args to config
    initial_config["gene_name"] = getattr(args, "gene_name", "") or ""
    initial_config["gene_file"] = getattr(args, "gene_file", None)

    # Merge VCF file path into config
    initial_config["vcf_file"] = vcf_file

    # Ensure all important config values are present
    # The config from cli.py has all merged values we need
    # Make sure we preserve all of them, not just a subset
    # Note: initial_config already has the full merged config if passed as dict

    # Only add args that might not be in config yet
    args_to_add = [
        "reference",
        "filters",
        "gzip_intermediates",
        "no_replacement",
        "late_filtering",
        "final_filter",  # Add final_filter to ensure it's transferred to config
        "threads",
        "keep_intermediates",
        "output_file",
        "no_links",
        "output_dir",
        "xlsx",
        "html_report",
        "inheritance_mode",  # Add inheritance_mode to config
        "calculate_inheritance",  # Add calculate_inheritance flag
        # Gene burden analysis arguments
        "case_samples_file",
        "control_samples_file",
        "case_samples",
        "control_samples",
        "case_phenotypes",
        "control_phenotypes",
        "case_phenotypes_file",
        "control_phenotypes_file",
        # IGV-related arguments
        "igv",
        "bam_mapping_file",
        "igv_reference",
        "igv_fasta",
        "igv_ideogram",
        "igv_flanking",
        "igv_max_allele_len_filename",
        "igv_hash_len_filename",
        "igv_max_variant_part_filename",
        # Memory management arguments
        "max_variants_per_gene",
    ]
    for key in args_to_add:
        if hasattr(args, key) and getattr(args, key) is not None and key not in initial_config:
            initial_config[key] = getattr(args, key)

    # Map igv argument to igv_enabled for consistency with old pipeline
    if hasattr(args, "igv"):
        initial_config["igv_enabled"] = args.igv

    # Create pipeline context
    context = PipelineContext(args=args, config=initial_config, workspace=workspace)

    # Initialize checkpoint system if enabled
    if initial_config.get("enable_checkpoint", False):
        from .checkpoint import PipelineState

        context.checkpoint_state = PipelineState(
            initial_config["output_dir"],
            enable_checksum=initial_config.get("checkpoint_checksum", False),
        )

        # Initialize checkpoint with current pipeline version and config
        pipeline_version = initial_config.get("pipeline_version", "refactored_pipeline")

        # Handle checkpoint state initialization/loading
        is_resuming = initial_config.get("resume", False) or initial_config.get("resume_from")
        if is_resuming:
            # Resume mode: try to load existing checkpoint state
            if context.checkpoint_state.load():
                logger.info("Loaded existing checkpoint state for resume")
            else:
                logger.error("Cannot resume: No checkpoint file found")
                raise ValueError("Resume requested but no checkpoint file found")
        else:
            # Normal mode: initialize fresh checkpoint state (overwrites existing)
            context.checkpoint_state.initialize(initial_config, pipeline_version)
            logger.info("Initialized fresh checkpoint state")

    # Build stages based on configuration
    stages = build_pipeline_stages(args)

    logger.info(f"Pipeline configured with {len(stages)} stages")

    # Create runner
    enable_checkpoints = getattr(args, "enable_checkpoint", False)
    max_workers = getattr(args, "threads", None)
    # Use process executor for better CPU-bound performance when multiple threads requested
    executor_type = "process" if max_workers and max_workers > 1 else "thread"
    runner = PipelineRunner(
        enable_checkpoints=enable_checkpoints,
        max_workers=max_workers,
        executor_type=executor_type,
        enable_stage_batching=True,
    )

    # Show execution plan in debug mode
    if log_level == "DEBUG":
        plan = runner.dry_run(stages)
        logger.debug("Execution plan:")
        for i, level in enumerate(plan):
            logger.debug(f"  Level {i}: {level}")

    try:
        # Run pipeline
        runner.run(stages, context)

        logger.info("Pipeline completed successfully!")

        # Clean up intermediates if requested
        if not getattr(args, "keep_intermediates", False):
            workspace.cleanup(keep_intermediates=False)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise
    finally:
        # Always clean temp files
        workspace.cleanup(keep_intermediates=True)


def main():
    """Show example main function demonstrating how to use the refactored pipeline."""
    # This would normally come from CLI parsing
    args = argparse.Namespace(
        vcf_file="example.vcf",
        gene_name="BRCA1 BRCA2",
        output_dir="output",
        threads=4,
        filter="QUAL >= 30",
        extract=["CHROM", "POS", "REF", "ALT", "QUAL"],
        replace_genotypes=True,
        html_report=True,
        xlsx=True,
        enable_checkpoint=True,
    )

    run_refactored_pipeline(args)


if __name__ == "__main__":
    main()

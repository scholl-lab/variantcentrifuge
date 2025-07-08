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
from .pipeline import compute_base_name
from .pipeline_core.context import PipelineContext
from .pipeline_core.runner import PipelineRunner
from .pipeline_core.workspace import Workspace
from .gene_bed import normalize_genes
from .stages.setup_stages import (
    ConfigurationLoadingStage,
    PhenotypeLoadingStage,
    ScoringConfigLoadingStage,
    PedigreeLoadingStage,
    AnnotationConfigLoadingStage,
    SampleConfigLoadingStage,
)
from .stages.processing_stages import (
    GeneBedCreationStage,
    VariantExtractionStage,
    ParallelVariantExtractionStage,
    MultiAllelicSplitStage,
    SnpSiftFilterStage,
    FieldExtractionStage,
    GenotypeReplacementStage,
    PhenotypeIntegrationStage,
    ExtraColumnRemovalStage,
)
from .stages.analysis_stages import (
    DataFrameLoadingStage,
    CustomAnnotationStage,
    InheritanceAnalysisStage,
    VariantScoringStage,
    StatisticsGenerationStage,
    GeneBurdenAnalysisStage,
    VariantAnalysisStage,
)
from .stages.output_stages import (
    VariantIdentifierStage,
    FinalFilteringStage,
    PseudonymizationStage,
    TSVOutputStage,
    ExcelReportStage,
    HTMLReportStage,
    IGVReportStage,
    MetadataGenerationStage,
    ArchiveCreationStage,
    ParallelReportGenerationStage,
)

logger = logging.getLogger(__name__)


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

    # Setup stages (can run in parallel)
    if args.phenotype_file:
        stages.append(PhenotypeLoadingStage())

    if hasattr(args, "scoring_config_path") and args.scoring_config_path:
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

    # Processing stages
    stages.append(GeneBedCreationStage())

    # Choose parallel or single-threaded extraction
    threads = getattr(args, "threads", 1)
    if threads > 1:
        stages.append(ParallelVariantExtractionStage())
    else:
        stages.append(VariantExtractionStage())

    # Optional processing stages
    if hasattr(args, "snpeff_split_by_transcript") and args.snpeff_split_by_transcript:
        stages.append(MultiAllelicSplitStage())

    # Add SnpSift filtering stage (it will check config for filters at runtime)
    if not getattr(args, "late_filtering", False):
        stages.append(SnpSiftFilterStage())

    # Field extraction
    stages.append(FieldExtractionStage())

    # Optional transformations
    # Default behavior: do genotype replacement unless --no-replacement is specified
    if not getattr(args, "no_replacement", False):
        stages.append(GenotypeReplacementStage())

    if args.phenotype_file:
        stages.append(PhenotypeIntegrationStage())

    if hasattr(args, "extra_columns_to_remove") and args.extra_columns_to_remove:
        stages.append(ExtraColumnRemovalStage())

    # Analysis stages
    stages.append(DataFrameLoadingStage())

    if any(
        [
            hasattr(args, "annotate_bed") and args.annotate_bed,
            hasattr(args, "annotate_gene_list") and args.annotate_gene_list,
            hasattr(args, "annotate_json_genes") and args.annotate_json_genes,
        ]
    ):
        stages.append(CustomAnnotationStage())

    if hasattr(args, "ped_file") and args.ped_file:
        stages.append(InheritanceAnalysisStage())
    elif hasattr(args, "inheritance_mode") and args.inheritance_mode:
        # Can do inheritance without PED (all samples as affected singletons)
        stages.append(InheritanceAnalysisStage())

    if hasattr(args, "scoring_config_path") and args.scoring_config_path:
        stages.append(VariantScoringStage())

    # Always run variant analysis to add standard columns
    stages.append(VariantAnalysisStage())

    # Add VAR_ID after variant analysis since it creates a new DataFrame
    stages.append(VariantIdentifierStage())

    if not getattr(args, "no_stats", False):
        stages.append(StatisticsGenerationStage())

    if hasattr(args, "perform_gene_burden") and args.perform_gene_burden:
        stages.append(GeneBurdenAnalysisStage())

    if getattr(args, "late_filtering", False) or getattr(args, "final_filter", None):
        stages.append(FinalFilteringStage())

    if hasattr(args, "pseudonymize") and args.pseudonymize:
        stages.append(PseudonymizationStage())

    stages.append(TSVOutputStage())

    # Report generation - use parallel stage if multiple reports
    report_count = sum(
        [
            getattr(args, "xlsx", False) or getattr(args, "excel", False),
            getattr(args, "html_report", False),
            getattr(args, "igv_report", False),
            not getattr(args, "no_metadata", False),
        ]
    )

    if report_count > 1:
        stages.append(ParallelReportGenerationStage())
    else:
        # Add individual report stages
        if getattr(args, "xlsx", False) or getattr(args, "excel", False):
            stages.append(ExcelReportStage())

        if getattr(args, "html_report", False):
            stages.append(HTMLReportStage())

        if getattr(args, "igv_report", False):
            stages.append(IGVReportStage())

        if not getattr(args, "no_metadata", False):
            stages.append(MetadataGenerationStage())

    if hasattr(args, "archive_results") and args.archive_results:
        stages.append(ArchiveCreationStage())

    return stages


def run_refactored_pipeline(args: argparse.Namespace) -> None:
    """Run the refactored pipeline with the new architecture.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments
    """
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
        f"Checking args.config: hasattr={hasattr(args, 'config')}, type={type(getattr(args, 'config', None)) if hasattr(args, 'config') else 'N/A'}"
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
    if "fields_to_extract" in initial_config:
        logger.debug(f"fields_to_extract value: {initial_config['fields_to_extract'][:100]}...")
    if "extract" in initial_config:
        logger.debug(f"extract fields count: {len(initial_config['extract'])}")

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
        "threads",
        "keep_intermediates",
        "output_file",
        "no_links",
        "output_dir",
        "xlsx",
        "html_report",
    ]
    for key in args_to_add:
        if hasattr(args, key) and getattr(args, key) is not None and key not in initial_config:
            initial_config[key] = getattr(args, key)

    # Create pipeline context
    context = PipelineContext(args=args, config=initial_config, workspace=workspace)

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

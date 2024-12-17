# variantcentrifuge/cli.py

# This script serves as the command-line interface for variantcentrifuge.
# It parses user arguments, loads configuration, validates inputs, and
# delegates execution to the main pipeline, producing filtered and analyzed
# variant results.

"""
CLI module for variantcentrifuge.

This module:
- Parses command line arguments.
- Sets up logging and configuration.
- Validates mandatory fields (reference, filters, and fields) before
  proceeding with analysis.
- Validates input files (VCF, phenotype, gene files).
- Invokes the main pipeline to perform variant filtering, extraction,
  optional genotype replacement, phenotype integration, analyses, and
  result production.
"""

import argparse
import sys
import logging
import datetime
import os

from . import __version__
from .config import load_config
from .validators import validate_mandatory_parameters, validate_vcf_file, validate_phenotype_file
from .pipeline import run_pipeline

logger = logging.getLogger("variantcentrifuge")


def main():
    """
    Main entry point for variantcentrifuge CLI.

    Steps:
    - Parse arguments.
    - Configure logging and load config.
    - Validate mandatory parameters (reference, filters, fields).
    - Validate input files (VCF, phenotype).
    - Update configuration with CLI parameters.
    - Run the pipeline.
    """
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s",
                        stream=sys.stderr)

    start_time = datetime.datetime.now()
    logger.info(f"Run started at {start_time.isoformat()}")

    parser = argparse.ArgumentParser(
        description="variantcentrifuge: Filter and process VCF files."
    )
    parser.add_argument("--version",
                        action="version",
                        version=f"variantcentrifuge {__version__}",
                        help="Show the current version and exit")
    parser.add_argument("--log-level",
                        choices=["DEBUG", "INFO", "WARN", "ERROR"],
                        default="INFO", help="Set the logging level")
    parser.add_argument("-c", "--config",
                        help="Path to configuration file", default=None)
    parser.add_argument("-g", "--gene-name",
                        help="Gene or list of genes of interest")
    parser.add_argument("-G", "--gene-file",
                        help="File containing gene names, one per line")
    parser.add_argument("-v", "--vcf-file", help="Input VCF file path",
                        required=True)
    parser.add_argument("-r", "--reference",
                        help="Reference database for snpEff", default=None)
    parser.add_argument("-f", "--filters",
                        help="Filters to apply in SnpSift filter")
    parser.add_argument("-e", "--fields",
                        help="Fields to extract with SnpSift extractFields")
    parser.add_argument("-o", "--output-file", nargs='?', const='stdout',
                        help="Final output file name or 'stdout'/'-' for stdout")
    parser.add_argument("--xlsx", action="store_true",
                        help="Convert final result to Excel")
    parser.add_argument("--no-replacement", action="store_true",
                        help="Skip genotype replacement step")
    parser.add_argument("--stats-output-file",
                        help="File to write analysis statistics")
    parser.add_argument("--perform-gene-burden", action="store_true",
                        help="Perform gene burden analysis")
    parser.add_argument("--output-dir", help="Directory to store intermediate "
                        "and final output files", default="output")
    parser.add_argument("--keep-intermediates", action="store_true",
                        default=True, help="Keep intermediate files.")
    parser.add_argument("--phenotype-file",
                        help="Path to phenotype file (.csv or .tsv)")
    parser.add_argument("--phenotype-sample-column",
                        help="Name of the column containing sample IDs in "
                             "phenotype file")
    parser.add_argument("--phenotype-value-column",
                        help="Name of the column containing phenotype values "
                             "in phenotype file")
    parser.add_argument("--no-stats", action="store_true", default=False,
                        help="Skip the statistics computation step.")

    parser.add_argument("--case-phenotypes",
                        help="Comma-separated HPO terms defining case group")
    parser.add_argument("--control-phenotypes",
                        help="Comma-separated HPO terms defining control group")
    parser.add_argument("--case-phenotypes-file",
                        help="File with HPO terms for case group")
    parser.add_argument("--control-phenotypes-file",
                        help="File with HPO terms for control group")

    parser.add_argument("--case-samples",
                        help="Comma-separated sample IDs defining the case group")
    parser.add_argument("--control-samples",
                        help="Comma-separated sample IDs defining the control group")
    parser.add_argument("--case-samples-file",
                        help="File with sample IDs for case group")
    parser.add_argument("--control-samples-file",
                        help="File with sample IDs for control group")

    parser.add_argument("--gene-burden-mode", choices=["samples", "alleles"],
                        default="alleles",
                        help="Mode for gene burden calculation: 'samples' or 'alleles'")
    parser.add_argument("--correction-method", choices=["fdr", "bonferroni"],
                        default="fdr",
                        help="Multiple testing correction method for gene burden test")

    args = parser.parse_args()

    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "ERROR": logging.ERROR
    }
    logging.getLogger("variantcentrifuge").setLevel(
        log_level_map[args.log_level]
    )

    logger.debug(f"CLI arguments: {args}")

    # Load configuration
    cfg = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")

    # Merge CLI params with config if needed
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")

    # Validate mandatory parameters
    validate_mandatory_parameters(reference, filters, fields)

    # Validate VCF file
    validate_vcf_file(args.vcf_file, logger)

    # Validate phenotype file if provided
    if args.phenotype_file and (args.phenotype_sample_column is None or
                                args.phenotype_value_column is None):
        logger.error(
            "If a phenotype file is provided, both --phenotype-sample-column "
            "and --phenotype-value-column must be specified."
        )
        sys.exit(1)
    if args.phenotype_file:
        validate_phenotype_file(
            args.phenotype_file,
            args.phenotype_sample_column,
            args.phenotype_value_column,
            logger
        )

    # Update cfg with CLI parameters
    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats
    cfg["gene_burden_mode"] = args.gene_burden_mode
    cfg["correction_method"] = args.correction_method

    # Hand off execution to pipeline
    run_pipeline(args, cfg, start_time)


if __name__ == "__main__":
    main()

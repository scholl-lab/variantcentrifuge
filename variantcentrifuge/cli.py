# File: variantcentrifuge/cli.py

import argparse
import sys
import logging
import datetime
import os
from typing import Optional, Dict, Any

from .version import __version__
from .config import load_config
from .validators import validate_mandatory_parameters, validate_vcf_file, validate_phenotype_file
from .pipeline import run_pipeline

logger = logging.getLogger("variantcentrifuge")


def main() -> None:
    """
    Main entry point for variantcentrifuge CLI.

    Steps:
        1. Parse arguments.
        2. Configure logging and load config.
        3. Validate mandatory parameters (reference, filters, fields).
        4. Validate input files (VCF, phenotype).
        5. Update configuration with CLI parameters.
        6. Run the pipeline.

    Changes for presets:
        - Added a --preset argument which allows specifying one or more predefined filters.
        - If multiple presets are chosen, they are combined with AND (&).
        - If user-specified filters are also provided, they are combined with these presets using AND.
        - Presets must be defined in config.json under "presets".

    Changes for sample substring removal:
        - Added a --remove-sample-substring argument which, if provided, removes the specified substring
          from all sample names extracted from the VCF before any comparisons or mappings.

    Changes for links:
        - Added a --no-links argument to disable adding link columns. By default, links are added.

    Changes for custom columns:
        - Added a --add-column argument which can be given multiple times to append columns with given
          headers (but blank values) to the final output.

    Changes for transcript filtering:
        - Added a --transcript-list and --transcript-file argument which allows specifying one or more
          transcript IDs (by comma-separated list or a file with one transcript ID per line).
        - If provided, these transcripts are used to construct a SnpSift filter expression that
          filters variants by EFF[*].TRID field.

    New changes for extra sample fields in genotype:
        - Added --append-extra-sample-fields, --extra-sample-fields, and
          --extra-sample-field-delimiter arguments to allow appending fields
          like DP, AD, etc. alongside each genotype.
    """
    # Initial basic logging setup to stderr before arguments are parsed,
    # so that we can log early messages like start time.
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr
    )

    start_time: datetime.datetime = datetime.datetime.now()
    logger.info(f"Run started at {start_time.isoformat()}")

    parser = argparse.ArgumentParser(
        description="variantcentrifuge: Filter and process VCF files."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"variantcentrifuge {__version__}",
        help="Show the current version and exit",
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARN", "ERROR"],
        default="INFO",
        help="Set the logging level",
    )
    parser.add_argument(
        "--log-file",
        help="Path to a file to write logs to (in addition to stderr)."
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Path to configuration file",
        default=None,
    )
    parser.add_argument(
        "-g",
        "--gene-name",
        help="Gene or list of genes of interest"
    )
    parser.add_argument(
        "-G",
        "--gene-file",
        help="File containing gene names, one per line"
    )
    parser.add_argument(
        "-v",
        "--vcf-file",
        help="Input VCF file path",
        required=True
    )
    parser.add_argument(
        "-r",
        "--reference",
        help="Reference database for snpEff",
        default=None,
    )
    parser.add_argument(
        "-f",
        "--filters",
        help="Filters to apply in SnpSift filter"
    )
    # New argument for presets:
    parser.add_argument(
        "--preset",
        action="append",
        help=(
            "Apply predefined filtering presets defined in the config file. "
            "Specify multiple times for multiple presets. They are combined with AND. "
            "If custom filters are also given, they are combined with these presets using AND."
        )
    )
    parser.add_argument(
        "-e",
        "--fields",
        help="Fields to extract with SnpSift extractFields"
    )
    parser.add_argument(
        "-o",
        "--output-file",
        nargs='?',
        const='stdout',
        help="Final output file name or 'stdout'/'-' for stdout"
    )
    parser.add_argument(
        "--xlsx",
        action="store_true",
        help="Convert final result to Excel"
    )
    parser.add_argument(
        "--no-replacement",
        action="store_true",
        help="Skip genotype replacement step"
    )
    parser.add_argument(
        "--stats-output-file",
        help="File to write analysis statistics"
    )
    parser.add_argument(
        "--perform-gene-burden",
        action="store_true",
        help="Perform gene burden analysis"
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to store intermediate and final output files",
        default="output"
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        default=True,
        help="Keep intermediate files."
    )
    parser.add_argument(
        "--phenotype-file",
        help="Path to phenotype file (.csv or .tsv)"
    )
    parser.add_argument(
        "--phenotype-sample-column",
        help="Name of the column containing sample IDs in phenotype file"
    )
    parser.add_argument(
        "--phenotype-value-column",
        help="Name of the column containing phenotype values in phenotype file"
    )
    parser.add_argument(
        "--no-stats",
        action="store_true",
        default=False,
        help="Skip the statistics computation step."
    )
    parser.add_argument(
        "--case-phenotypes",
        help="Comma-separated HPO terms defining case group"
    )
    parser.add_argument(
        "--control-phenotypes",
        help="Comma-separated HPO terms defining control group"
    )
    parser.add_argument(
        "--case-phenotypes-file",
        help="File with HPO terms for case group"
    )
    parser.add_argument(
        "--control-phenotypes-file",
        help="File with HPO terms for control group"
    )
    parser.add_argument(
        "--case-samples",
        help="Comma-separated sample IDs defining the case group"
    )
    parser.add_argument(
        "--control-samples",
        help="Comma-separated sample IDs defining the control group"
    )
    parser.add_argument(
        "--case-samples-file",
        help="File with sample IDs for case group"
    )
    parser.add_argument(
        "--control-samples-file",
        help="File with sample IDs for control group"
    )
    parser.add_argument(
        "--gene-burden-mode",
        choices=["samples", "alleles"],
        default="alleles",
        help="Mode for gene burden calculation: 'samples' or 'alleles'"
    )
    parser.add_argument(
        "--correction-method",
        choices=["fdr", "bonferroni"],
        default="fdr",
        help="Multiple testing correction method for gene burden test"
    )
    parser.add_argument(
        "--html-report",
        action="store_true",
        help="Generate an interactive HTML report with sortable variant tables and summary plots."
    )
    parser.add_argument(
        "--igv",
        action="store_true",
        help="Enable IGV.js integration for genomic visualization."
    )
    parser.add_argument(
        "--bam-mapping-file",
        help="Path to a TSV or CSV file mapping sample IDs to BAM files (sample_id,bam_path)."
    )
    parser.add_argument(
        "--igv-reference",
        help="Genome reference identifier for IGV (e.g., 'hg19' or 'hg38'). Required if --igv is enabled."
    )

    # New argument for sample name substring removal
    parser.add_argument(
        "--remove-sample-substring",
        help="If provided, this substring will be removed from all sample names found in the VCF."
    )

    # Added argument to control adding links
    parser.add_argument(
        "--no-links",
        action="store_true",
        default=False,
        help="Disable adding link columns to the final output (links are added by default)."
    )

    # >>> New threads argument
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads to use for bcftools and related operations."
    )

    # >>> New argument for adding named columns
    parser.add_argument(
        "--add-column",
        action="append",
        default=[],
        help="Append a named blank column to the final output. Repeat for multiple columns."
    )

    # >>> New argument to split SNPeff lines
    parser.add_argument(
        "--split-snpeff-lines",
        action="store_true",
        default=False,
        help=(
            "If set, run the SNPeff annotation splitter step on the extracted variants VCF "
            "before applying SNPsift filters. This splits multiple EFF/ANN annotations into "
            "separate VCF lines."
        )
    )

    # >>> New arguments for transcript filtering
    parser.add_argument(
        "--transcript-list",
        help=(
            "Comma-separated list of transcript IDs to filter for. "
            "For example: NM_007294.4,NM_000059.4"
        )
    )
    parser.add_argument(
        "--transcript-file",
        help=(
            "Path to a file containing transcript IDs, one per line. "
            "For example: a file where each line is a transcript like NM_007294.4"
        )
    )

    # >>> NEW arguments for genotype filtering
    parser.add_argument(
        "--genotype-filter",
        help=(
            "Specify genotype filter(s), e.g. 'het', 'hom', 'comp_het' or "
            "comma-separated combination (e.g. 'het,comp_het'). Only samples/lines "
            "fulfilling these genotypes are kept (unless per-gene file overrides)."
        )
    )
    parser.add_argument(
        "--gene-genotype-file",
        help=(
            "Path to a file specifying per-gene genotype rules. It must contain at least two columns: "
            "'GENE' and 'GENOTYPES' (a comma-separated list of genotype filters, e.g. het,comp_het). "
            "If a gene is listed here, it overrides any global --genotype-filter."
        )
    )

    # >>> NEW arguments for appending extra sample fields
    parser.add_argument(
        "--append-extra-sample-fields",
        action="store_true",
        default=False,
        help=(
            "If set, append extra sample fields (e.g. DP, AD) to the genotype in the final output, "
            "e.g. sample(0/1:DP_value:AD_value). Requires that these fields are included in "
            "--fields or cfg['fields_to_extract']. Defaults to False."
        ),
    )
    parser.add_argument(
        "--extra-sample-fields",
        default="",
        help=(
            "Comma-separated list of additional sample-level fields (e.g. DP,AD) to append next "
            "to each genotype in the final output if --append-extra-sample-fields is set. "
            "These fields must also be present in the extracted TSV."
        ),
    )
    parser.add_argument(
        "--extra-sample-field-delimiter",
        default=":",
        help=(
            "Delimiter to separate genotype from extra sample fields (and from each other). "
            "Default is ':'."
        ),
    )

    args: argparse.Namespace = parser.parse_args()

    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    logging.getLogger("variantcentrifuge").setLevel(log_level_map[args.log_level])

    # If a log file is specified, add a file handler in addition to stderr
    if args.log_file:
        fh = logging.FileHandler(args.log_file)
        fh.setLevel(log_level_map[args.log_level])
        fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(fh)
        logger.debug(f"Logging to file enabled: {args.log_file}")

    logger.debug(f"CLI arguments: {args}")

    # Load configuration
    cfg: Dict[str, Any] = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")

    reference: Optional[str] = args.reference or cfg.get("reference")
    filters: Optional[str] = args.filters or cfg.get("filters")
    fields: Optional[str] = args.fields or cfg.get("fields_to_extract")

    # Handle presets if provided
    preset_dict = cfg.get("presets", {})
    if args.preset:
        chosen_presets = []
        for p in args.preset:
            p_str = p.strip()
            if p_str in preset_dict:
                chosen_presets.append(f"({preset_dict[p_str]})")
            else:
                logger.error(f"Preset '{p_str}' not found in config.")
                sys.exit(1)
        if chosen_presets:
            combined_preset_filter = " & ".join(chosen_presets)
            if filters:
                filters = f"({combined_preset_filter}) & ({filters})"
            else:
                filters = combined_preset_filter

    # Validate mandatory parameters
    validate_mandatory_parameters(reference, filters, fields)
    validate_vcf_file(args.vcf_file, logger)

    if args.phenotype_file and (
        args.phenotype_sample_column is None or args.phenotype_value_column is None
    ):
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

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats
    cfg["gene_burden_mode"] = args.gene_burden_mode
    cfg["correction_method"] = args.correction_method

    # IGV parameters
    cfg["igv_enabled"] = args.igv
    cfg["bam_mapping_file"] = args.bam_mapping_file
    cfg["igv_reference"] = args.igv_reference

    # If IGV is enabled, check that bam_mapping_file and igv_reference are provided
    if args.igv and (not args.bam_mapping_file or not args.igv_reference):
        logger.error("For IGV integration, --bam-mapping-file and --igv-reference must be provided.")
        sys.exit(1)

    # Update filters and fields in cfg after handling presets
    cfg["filters"] = filters
    cfg["fields_to_extract"] = fields

    # Store remove_sample_substring parameter if given
    if args.remove_sample_substring:
        cfg["remove_sample_substring"] = args.remove_sample_substring
    else:
        cfg["remove_sample_substring"] = None

    # Store no_links parameter
    cfg["no_links"] = args.no_links

    # >>> Store threads in cfg
    cfg["threads"] = args.threads

    # >>> Store the transcript arguments in cfg so pipeline can see them
    cfg["transcript_list"] = args.transcript_list  # e.g. "NM_007294.4,NM_000059.4"
    cfg["transcript_file"] = args.transcript_file  # path to file with transcripts

    # >>> Store genotype-filter arguments
    # (already in 'args', no further validation required here)
    cfg["genotype_filter"] = args.genotype_filter
    cfg["gene_genotype_file"] = args.gene_genotype_file

    # >>> Store new arguments for extra sample fields
    cfg["append_extra_sample_fields"] = args.append_extra_sample_fields
    # Convert CSV string to a list
    cfg["extra_sample_fields"] = [
        f.strip() for f in args.extra_sample_fields.split(",") if f.strip()
    ]
    cfg["extra_sample_field_delimiter"] = args.extra_sample_field_delimiter

    # Finally, run the pipeline
    run_pipeline(args, cfg, start_time)

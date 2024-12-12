# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

"""
Command-line interface (CLI) module.

This module defines the main entry point for the variantcentrifuge CLI,
handling argument parsing, orchestrating steps with temporary files, and
following logic similar to the original bash script.
"""

import argparse
import sys
import tempfile
from . import __version__
from .config import load_config
from .gene_bed import get_gene_bed
from .filters import extract_variants, apply_snpsift_filter
from .extractor import extract_fields
from .replacer import replace_genotypes
from .phenotype_filter import filter_phenotypes
from .converter import convert_to_excel
from .analyze_variants import analyze_variants
from .utils import log_message, check_external_tools, set_log_level


def main():
    parser = argparse.ArgumentParser(
        description="variantcentrifuge: Filter and process VCF files."
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"variantcentrifuge {__version__}",
        help="Show the current version and exit"
    )
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARN", "ERROR"],
        default="INFO",
        help="Set the logging level (DEBUG, INFO, WARN, ERROR)"
    )
    parser.add_argument(
        "-c", "--config", help="Path to configuration file", default=None
    )
    parser.add_argument(
        "-g", "--gene-name", help="Gene or list of genes of interest",
        required=True
    )
    parser.add_argument(
        "-v", "--vcf-file", help="Input VCF file path", required=True
    )
    parser.add_argument(
        "-r", "--reference", help="Reference database for snpEff",
        default=None
    )
    parser.add_argument(
        "-f", "--filters", help="Filters to apply in SnpSift filter"
    )
    parser.add_argument(
        "-e", "--fields", help="Fields to extract with SnpSift extractFields"
    )
    parser.add_argument(
        "-s", "--samples-file", help="Samples file for genotype replacement"
    )
    parser.add_argument(
        "-o", "--output-file", help="Output file name"
    )
    parser.add_argument(
        "--xlsx", action="store_true", help="Convert final result to Excel"
    )
    parser.add_argument(
        "--no-replacement", action="store_true",
        help="Skip genotype replacement step"
    )
    parser.add_argument(
        "--stats-output-file", help="File to write analysis statistics"
    )
    parser.add_argument(
        "--perform-gene-burden", action="store_true",
        help="Perform gene burden analysis"
    )

    args = parser.parse_args()
    set_log_level(args.log_level)
    log_message("DEBUG", f"CLI arguments: {args}")

    # Check external tools
    check_external_tools()

    # Load config
    cfg = load_config(args.config)
    log_message("DEBUG", f"Configuration loaded: {cfg}")
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")
    output_file = args.output_file or cfg.get("output_file", "variants.tsv")

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["stats_output_file"] = args.stats_output_file

    # Generate BED file from gene
    log_message("DEBUG", "Starting gene BED extraction.")
    bed_file = get_gene_bed(
        reference,
        args.gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True)
    )
    log_message("DEBUG", f"Gene BED created at: {bed_file}")

    # Extract variants with bcftools
    log_message("DEBUG", "Extracting variants with bcftools...")
    variant_file = extract_variants(args.vcf_file, bed_file, cfg)
    log_message("DEBUG", f"Variants extracted to {variant_file}. Applying filters...")

    # Apply SnpSift filter
    filtered_file = apply_snpsift_filter(variant_file, filters, cfg)
    log_message("DEBUG", f"Filter applied. Extracting fields to TSV...")

    # Extract fields
    extracted_file = extract_fields(filtered_file, fields, cfg)
    log_message("DEBUG", f"Fields extracted to {extracted_file}")

    # Genotype replacement if not skipped
    if not args.no_replacement:
        log_message("DEBUG", "Replacing genotypes if sample info is available...")
        replaced_data_file = tempfile.mktemp(suffix=".tsv")
        # replace_genotypes returns an iterator of lines, so let's just read/write directly
        # Adjust replace_genotypes to produce a file instead of lines:
        with open(extracted_file, "r", encoding="utf-8") as inp, open(replaced_data_file, "w", encoding="utf-8") as out:
            for line in replace_genotypes(inp, args.samples_file, cfg):
                out.write(line + "\n")
        log_message("DEBUG", f"Genotype replacement done, results in {replaced_data_file}")
    else:
        replaced_data_file = extracted_file
        log_message("DEBUG", "Genotype replacement skipped.")

    # Phenotype filtering if enabled
    if cfg.get("use_phenotype_filtering", False):
        log_message("DEBUG", "Applying phenotype filtering...")
        phenotype_file = tempfile.mktemp(suffix=".tsv")
        with open(replaced_data_file, "r", encoding="utf-8") as inp, open(phenotype_file, "w", encoding="utf-8") as out:
            for line in filter_phenotypes(inp, cfg):
                out.write(line + "\n")
        replaced_data_file = phenotype_file
        log_message("DEBUG", "Phenotype filtering complete.")

    # Analyze variants if requested
    log_message("DEBUG", "Analyzing variants if requested...")
    analyzed_data_file = tempfile.mktemp(suffix=".tsv")
    # analyze_variants returns an iterator of lines as well
    with open(replaced_data_file, "r", encoding="utf-8") as inp, open(analyzed_data_file, "w", encoding="utf-8") as out:
        for line in analyze_variants(inp, cfg):
            out.write(line + "\n")
    log_message("DEBUG", f"Variant analysis complete, results in {analyzed_data_file}")

    # Write final output
    final_output = output_file if output_file else "/dev/stdout"
    log_message("DEBUG", f"Writing final output to {final_output}")
    with open(analyzed_data_file, "r", encoding="utf-8") as inp, open(final_output, "w", encoding="utf-8") as out:
        for line in inp:
            out.write(line + "\n")

    if args.xlsx:
        log_message("DEBUG", "Converting final output to Excel...")
        convert_to_excel(final_output, cfg)
        log_message("DEBUG", "Excel conversion complete.")

    log_message("INFO", f"Processing completed. Output saved to {final_output}")

if __name__ == "__main__":
    main()

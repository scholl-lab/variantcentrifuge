# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

"""
Command-line interface (CLI) module.

This module defines the main entry point for the variantcentrifuge CLI,
handling argument parsing and orchestrating the pipeline steps.
"""

import argparse
import sys
from . import __version__
from .config import load_config
from .gene_bed import get_gene_bed
from .filters import extract_variants, apply_snpsift_filter
from .extractor import extract_fields
from .replacer import replace_genotypes
from .phenotype_filter import filter_phenotypes
from .converter import convert_to_excel
from .analyze_variants import analyze_variants
from .utils import log_message, check_external_tools

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

    # Check external tool availability before proceeding
    check_external_tools()

    # Load configuration and override defaults if provided
    cfg = load_config(args.config)
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")
    output_file = args.output_file or cfg.get("output_file", "variants.tsv")

    # Update cfg with runtime arguments
    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["stats_output_file"] = args.stats_output_file

    # Process gene BED
    bed_file = get_gene_bed(
        reference,
        args.gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True)
    )

    # Extract variants
    variant_stream = extract_variants(args.vcf_file, bed_file, cfg)

    # Apply SnpSift filters
    filtered_stream = apply_snpsift_filter(variant_stream, filters, cfg)

    # Extract fields
    extracted_data = extract_fields(filtered_stream, fields, cfg)

    # Replace genotypes if not skipped
    if not args.no_replacement:
        replaced_data = replace_genotypes(extracted_data, args.samples_file, cfg)
    else:
        replaced_data = extracted_data

    # Apply phenotype filtering if configured
    if cfg.get("use_phenotype_filtering", False):
        replaced_data = filter_phenotypes(replaced_data, cfg)

    # Analyze variants if requested
    analyzed_data = analyze_variants(replaced_data, cfg)

    # Write output
    with open(output_file, "w", encoding="utf-8") as out_f:
        for line in analyzed_data:
            out_f.write(line + "\n")

    # Convert to Excel if requested
    if args.xlsx:
        convert_to_excel(output_file, cfg)

    log_message("INFO", f"Processing completed. Output saved to {output_file}")

if __name__ == "__main__":
    main()

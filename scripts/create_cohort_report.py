#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
create_cohort_report.py

This script aggregates results from multiple single-sample VariantCentrifuge reports
into one comprehensive, interactive cohort-level report.

Author: VariantCentrifuge Team
"""

import os
import re
import sys
import json
import glob
import argparse
import logging
import pandas as pd
from jinja2 import Environment, FileSystemLoader

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("create_cohort_report")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create a cohort report from multiple VariantCentrifuge sample reports."
    )

    parser.add_argument(
        "--input-pattern",
        required=True,
        help='Glob pattern to find input TSV files (e.g., "/path/to/samples/*/variants.tsv")',
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory where the cohort report will be saved"
    )
    parser.add_argument(
        "--report-name",
        default="Variant Cohort Analysis",
        help='Name for the cohort report (default: "Variant Cohort Analysis")',
    )
    parser.add_argument(
        "--sample-regex",
        default=r".*?([^/\\]+)(?:/|\\)variants\.tsv$",
        help="Regex pattern to extract sample ID from input file paths",
    )

    return parser.parse_args()


def find_input_files(input_pattern):
    """Find all input TSV files using the provided glob pattern."""
    files = glob.glob(input_pattern, recursive=True)

    if not files:
        logger.error(f"No input files found matching pattern: {input_pattern}")
        sys.exit(1)

    logger.info(f"Found {len(files)} input files")
    return files


def extract_sample_id(file_path, regex_pattern):
    """Extract sample ID from file path using the provided regex pattern."""
    match = re.search(regex_pattern, file_path)
    if match and match.group(1):
        return match.group(1)
    # Fallback: use filename without extension if regex doesn't match
    return os.path.splitext(os.path.basename(file_path))[0]


def aggregate_data(input_files, sample_regex):
    """Aggregate data from multiple TSV files into a single DataFrame."""
    logger.info("Aggregating data from input files...")

    dfs = []

    for file_path in input_files:
        try:
            # Read the TSV file
            df = pd.read_csv(file_path, sep="\t")

            # Extract sample ID and add as a column
            sample_id = extract_sample_id(file_path, sample_regex)
            df["SampleID"] = sample_id

            # Look for sample-specific allele frequency column (may be named differently)
            af_col = None
            for col in df.columns:
                if "AF" in col and col != "AF_EXAC" and col != "AF_GNOMAD" and col != "AF_1000G":
                    af_col = col
                    break

            # Rename sample-specific allele frequency column if found
            if af_col and af_col != "SampleAF":
                df.rename(columns={af_col: "SampleAF"}, inplace=True)

            dfs.append(df)
            logger.debug(f"Processed {file_path} (sample: {sample_id})")

        except Exception as e:
            logger.warning(f"Error processing {file_path}: {str(e)}")

    if not dfs:
        logger.error("No valid data found in input files")
        sys.exit(1)

    # Concatenate all DataFrames
    master_df = pd.concat(dfs, ignore_index=True)
    logger.info(f"Aggregated {len(master_df)} total variants from {len(dfs)} samples")

    return master_df


def clean_data(df):
    """Clean and prepare the aggregated data."""
    logger.info("Cleaning and preparing data...")

    # Ensure all relevant numeric columns are properly typed
    numeric_cols = [
        "SampleAF",
        "QUAL",
        "DP",
        "AF_EXAC",
        "AF_GNOMAD",
        "AF_1000G",
        "IMPACT_SEVERITY",
        "CADD_PHRED",
    ]

    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Clean up any NaN values in string columns
    string_cols = ["Gene", "GeneID", "IMPACT", "Effect", "AAChange", "SampleID"]
    for col in string_cols:
        if col in df.columns:
            df[col] = df[col].fillna("Unknown")

    # Sort by gene and impact severity
    if "Gene" in df.columns:
        df = df.sort_values(by=["Gene", "SampleID"])

    return df


def compute_statistics(df):
    """Compute summary statistics for the cohort."""
    logger.info("Computing summary statistics...")

    # Global summary
    total_variants = len(df)
    unique_samples = df["SampleID"].nunique()
    unique_genes = df["Gene"].nunique() if "Gene" in df.columns else 0

    # Per-gene statistics
    gene_summary = None
    if "Gene" in df.columns:
        # Count variants per gene
        gene_variants = df.groupby("Gene").size().reset_index(name="VariantCount")

        # Count unique samples per gene
        gene_samples = df.groupby("Gene")["SampleID"].nunique().reset_index(name="SampleCount")

        # Merge the two statistics
        gene_summary = pd.merge(gene_variants, gene_samples, on="Gene")

        # Sort by number of samples and then by number of variants
        gene_summary = gene_summary.sort_values(by=["SampleCount", "VariantCount"], ascending=False)

    summary = {
        "total_variants": int(total_variants),
        "unique_samples": int(unique_samples),
        "unique_genes": int(unique_genes),
        "gene_summary": gene_summary.to_dict(orient="records") if gene_summary is not None else [],
    }

    return summary


def save_data(df, summary, output_dir):
    """Save the prepared data as JSON files."""
    logger.info("Saving data to JSON files...")

    # Create data directory if it doesn't exist
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    # Save variants data
    variants_file = os.path.join(data_dir, "variants.json")
    with open(variants_file, "w") as f:
        # Convert to JSON records with date_format='iso' for proper datetime handling
        variants_json = df.to_json(orient="records", date_format="iso")
        f.write(variants_json)

    # Save summary data
    summary_file = os.path.join(data_dir, "summary.json")
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Data saved to {data_dir}")


def create_report(output_dir, report_name, df, summary):
    """Create the HTML report using the template and data."""
    logger.info("Creating HTML report...")

    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set up Jinja2 environment
    env = Environment(loader=FileSystemLoader(os.path.join(script_dir, "templates")))
    template = env.get_template("cohort_report_template.html")

    # Convert DataFrame and summary to JSON strings for embedding in HTML
    variants_json = df.to_json(orient="records", date_format="iso")
    summary_json = json.dumps(summary)

    # Render the template with embedded JSON data
    html_content = template.render(
        report_name=report_name,
        generation_date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        variants_json=variants_json,
        summary_json=summary_json
    )

    # Write the HTML to a file
    report_file = os.path.join(output_dir, "cohort_report.html")
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    logger.info(f"HTML report created: {report_file}")
    return report_file


def main():
    """Main function to execute the cohort report creation."""
    # Parse command line arguments
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Find input files
    input_files = find_input_files(args.input_pattern)

    # Aggregate data from input files
    df = aggregate_data(input_files, args.sample_regex)

    # Clean and prepare data
    df = clean_data(df)

    # Compute statistics
    summary = compute_statistics(df)

    # Create HTML report with embedded data (no need to save external JSON files)
    report_file = create_report(args.output_dir, args.report_name, df, summary)

    logger.info("Cohort report creation complete!")
    logger.info(f"Report available at: {report_file}")


if __name__ == "__main__":
    main()

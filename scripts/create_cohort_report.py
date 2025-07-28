#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Script to create a cohort report from multiple sample TSV files.

create_cohort_report.py

This script aggregates results from multiple single-sample VariantCentrifuge reports
into one comprehensive, interactive cohort-level report.

Author: VariantCentrifuge Team
"""

import argparse
import glob
import json
import logging
import os
import re
import sys

import pandas as pd
from jinja2 import Environment, FileSystemLoader
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter

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
    logger.debug(
        f"Attempting to extract sample ID from: '{file_path}' using regex: '{regex_pattern}'"
    )
    match = re.search(regex_pattern, file_path)
    if match:
        # Ensure group 1 exists and is not None
        if len(match.groups()) > 0 and match.group(1):
            sample_id = match.group(1)
            logger.info(
                f"Regex match successful. Extracted SampleID: '{sample_id}' from '{file_path}'"
            )
            return sample_id
        else:
            # This case handles if regex matches but capture group 1 is missing or empty
            logger.warning(
                f"Regex matched for '{file_path}', but capture group 1 was missing or empty. "
                f"Match groups: {match.groups()}. Regex: '{regex_pattern}'"
            )
    else:
        logger.warning(
            f"Regex pattern '{regex_pattern}' did not match for file path: '{file_path}'."
        )

    # Fallback: use the name of the parent directory of the file_path
    # This is often more robust if 'variants.tsv' (or similar) is consistently named
    # and located within a sample-specific directory.
    try:
        parent_dir_name = os.path.basename(os.path.dirname(file_path))
        if parent_dir_name:  # Ensure parent_dir_name is not empty
            logger.info(
                f"Using fallback: parent directory name '{parent_dir_name}' as SampleID "
                f"for '{file_path}'"
            )
            return parent_dir_name
        else:  # Handle cases like file_path being in root directory, e.g. "/variants.tsv"
            logger.warning(
                f"Fallback to parent directory name failed for '{file_path}' "
                f"(parent_dir_name is empty). Using filename as last resort."
            )
    except Exception as e:
        logger.error(f"Error during fallback parent directory extraction for '{file_path}': {e}")

    # Last resort fallback: filename without extension
    final_fallback_id = os.path.splitext(os.path.basename(file_path))[0]
    logger.info(
        f"Using last resort fallback: filename without extension '{final_fallback_id}' "
        f"as SampleID for '{file_path}'"
    )
    return final_fallback_id


def aggregate_data(input_files, sample_regex):
    """Aggregate data from multiple TSV files into a single DataFrame."""
    logger.info(
        f"Starting data aggregation from {len(input_files)} input file(s) using pattern "
        f"for sample_regex: '{sample_regex}'"
    )

    dfs = []
    collected_sample_ids = []  # To track all extracted IDs

    for file_path in input_files:
        logger.info(f"Processing file: {file_path}")  # More prominent log for each file
        try:
            # Use low_memory=False to suppress dtype warnings for variant files
            df = pd.read_csv(file_path, sep="\t", low_memory=False)
            logger.debug(f"Successfully read {file_path}. Columns: {df.columns.tolist()}")

            sample_id = extract_sample_id(
                file_path, sample_regex
            )  # Calls the updated extract_sample_id
            df["SampleID"] = sample_id
            collected_sample_ids.append(sample_id)
            logger.info(f"Assigned SampleID '{sample_id}' to data from {file_path}")

            # Look for sample-specific allele frequency column (may be named differently)
            af_col = None
            for col in df.columns:
                # Ensure 'AF' is part of the column name, and it's not one of the
                # standard population AF columns
                if "AF" in col and col not in [
                    "AF_EXAC",
                    "AF_GNOMAD",
                    "AF_1000G",
                    "SampleAF",  # Added SampleAF to avoid re-processing
                ]:
                    af_col = col
                    break

            if af_col:  # If a suitable AF column is found
                if "SampleAF" not in df.columns:  # And SampleAF doesn't already exist
                    df.rename(columns={af_col: "SampleAF"}, inplace=True)
                    msg = (
                        f"Renamed sample-specific AF column '{af_col}' to 'SampleAF' "
                        f"for sample '{sample_id}' from '{file_path}'"
                    )
                    logger.debug(msg)
                elif (
                    af_col != "SampleAF"
                ):  # If SampleAF exists but we found another candidate (e.g. TUMOR_AF)
                    logger.debug(
                        f"Column 'SampleAF' already exists. Did not rename '{af_col}' "
                        f"for sample '{sample_id}' from '{file_path}'"
                    )

            # Ensure the Gene column exists (case-insensitive check and rename)
            if "Gene" not in df.columns:
                found_gene_col_original_case = None
                for col_idx, col_name in enumerate(df.columns):
                    if col_name.lower() == "gene":
                        found_gene_col_original_case = col_name
                        break
                if found_gene_col_original_case:
                    df.rename(columns={found_gene_col_original_case: "Gene"}, inplace=True)
                    logger.debug(
                        f"Renamed column '{found_gene_col_original_case}' to 'Gene' "
                        f"for sample '{sample_id}' from '{file_path}'"
                    )
                else:
                    logger.warning(
                        f"'Gene' column (or case-insensitive 'gene') not found in {file_path} "
                        f"for sample '{sample_id}'. Adding 'Unknown' as placeholder."
                    )
                    df["Gene"] = "Unknown"  # Add placeholder if no gene column

            dfs.append(df)

        except Exception as e:
            logger.error(
                f"Failed to process {file_path}. Error: {str(e)}", exc_info=True
            )  # exc_info=True for traceback

    if not dfs:
        logger.critical(
            "CRITICAL: No dataframes were created after processing all input files. "
            "This means no valid data could be read or processed. Exiting."
        )
        sys.exit(1)

    master_df = pd.concat(dfs, ignore_index=True)

    unique_extracted_samples = sorted(
        list(set(collected_sample_ids))
    )  # Get sorted list of unique IDs
    logger.info(f"Aggregation complete. Total variants aggregated: {len(master_df)}.")
    logger.info(
        f"All SampleIDs extracted during aggregation (includes duplicates if any): "
        f"{collected_sample_ids}"
    )
    logger.info(
        f"Number of unique SampleIDs successfully extracted: {len(unique_extracted_samples)}. "
        f"Unique IDs found: {unique_extracted_samples}"
    )

    if "SampleID" in master_df.columns:
        df_unique_samples_count = master_df["SampleID"].nunique()
        df_unique_samples_list = sorted(list(master_df["SampleID"].unique()))
        logger.info(
            f"Verification: Unique SampleIDs present in the final aggregated DataFrame: "
            f"{df_unique_samples_count}. IDs: {df_unique_samples_list}"
        )
        if (
            len(unique_extracted_samples) != df_unique_samples_count
            or unique_extracted_samples != df_unique_samples_list
        ):
            logger.warning(
                "Potential Mismatch Alert! The set of unique SampleIDs extracted during "
                "file processing "
                f"({len(unique_extracted_samples)}: {unique_extracted_samples}) "
                "differs from the unique SampleIDs found in the final DataFrame "
                f"({df_unique_samples_count}: {df_unique_samples_list}). "
                "This could indicate an issue in data concatenation or SampleID assignment."
            )
    else:
        logger.error(
            "CRITICAL: 'SampleID' column is MISSING from the final aggregated DataFrame. "
            "Sample identification will fail."
        )

    return master_df


def clean_data(df):
    """Clean and prepare the aggregated data."""
    logger.info("Cleaning and preparing data...")

    # Ensure SampleID is the first column
    if "SampleID" in df.columns:
        # Get all columns except SampleID
        other_columns = [col for col in df.columns if col != "SampleID"]
        # Reorder columns with SampleID first
        column_order = ["SampleID"] + other_columns
        df = df[column_order]
        logger.info("Reordered columns with SampleID as first column")

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


def detect_link_columns(df):
    """Detect columns that contain URLs or should be treated as links."""
    link_columns = {}

    if df.empty:
        return link_columns

    for col in df.columns:
        if col.lower() in ["url", "link", "href"]:
            link_columns[col] = {"type": "standard", "display_text": "Link"}
            continue

        # Check if column contains URLs
        sample_values = df[col].dropna().astype(str).head(10)
        url_count = sum(
            1 for val in sample_values if val.startswith(("http://", "https://", "ftp://"))
        )

        # If more than 50% of non-null values are URLs, treat as link column
        if len(sample_values) > 0 and url_count / len(sample_values) > 0.5:
            link_columns[col] = {"type": "standard", "display_text": col.replace("_", " ")}

    logger.info(f"Detected link columns: {list(link_columns.keys())}")
    return link_columns


def prepare_column_metadata(df):
    """Prepare column metadata for template rendering."""
    link_columns = detect_link_columns(df)

    # Use the same default hidden columns as individual reports from config.json
    default_hidden_columns = [
        "QUAL",
        "AC",
        "FEATUREID",
        "AA_POS",
        "AA_LEN",
        "NMD_PERC",
        "dbNSFP_ALFA_Total_AC",
        "dbNSFP_REVEL_score",
        "splice_dbscSNV_rf_score",
        "splice_dbscSNV_ada_score",
        "splice_spidex_dpsi_zscore",
        "dbNSFP_clinvar_clnsig",
        "dbNSFP_gnomAD_exomes_AC",
        "dbNSFP_gnomAD_genomes_AC",
        "hgmd_CLASS",
        "proband_count",
        "control_count",
        "proband_variant_count",
        "control_variant_count",
        "proband_allele_count",
        "control_allele_count",
        "proband_homozygous_count",
        "control_homozygous_count",
    ]

    # Use the same hover-expand columns as individual reports from config.json
    truncate_columns = [
        "VAR_ID",
        "HGVS_P",
        "HGVS_C",
        "ID",
        "REF",
        "ALT",
        "EFFECT",
        "IMPACT",
        "FEATUREID",
        "dbNSFP_clinvar_clnsig",
        "ClinVar_CLNSIG",
        "GT",
    ]

    column_metadata = []
    for col in df.columns:
        is_link = col in link_columns
        is_hidden = col in default_hidden_columns
        needs_truncation = col in truncate_columns

        metadata = {
            "name": col,
            "display_name": col.replace("_", " "),
            "is_link": is_link,
            "link_info": link_columns.get(col, {}),
            "is_hidden_default": is_hidden,
            "needs_truncation": needs_truncation,
            "max_width": 250 if col == "GT" else (120 if needs_truncation else None),
        }
        column_metadata.append(metadata)

    return column_metadata


def create_excel_report(df, output_dir):
    """Create an Excel report with clickable links."""
    logger.info("Creating Excel report...")

    # Create Excel file path
    excel_file = os.path.join(output_dir, "cohort_variants.xlsx")

    # Write DataFrame to Excel
    df.to_excel(excel_file, index=False, sheet_name="Cohort_Variants")

    # Add formatting and clickable links
    finalize_cohort_excel(excel_file, df)

    logger.info(f"Excel report created: {excel_file}")
    return excel_file


def finalize_cohort_excel(xlsx_file, df):
    """Apply formatting and create clickable links in the Excel file."""
    from openpyxl.styles import Font

    wb = load_workbook(xlsx_file)
    ws = wb.active

    # Freeze the top row
    ws.freeze_panes = "A2"

    # Enable auto-filter
    max_col_letter = get_column_letter(ws.max_column)
    ws.auto_filter.ref = f"A1:{max_col_letter}1"

    logger.info("Adding clickable links to Excel file...")

    # Define link templates (similar to individual reports)
    link_templates = {
        "SpliceAI": (
            "https://spliceailookup.broadinstitute.org/"
            "#variant={CHROM}-{POS}-{REF}-{ALT}&hg=19&bc=basic&distance=500&mask=0&ra=0"
        ),
        "Franklin": (
            "https://franklin.genoox.com/clinical-db/variant/snp/{CHROM}-{POS}-{REF}-{ALT}-hg19"
        ),
        "Varsome": "https://varsome.com/variant/hg19/{CHROM}-{POS}-{REF}-{ALT}",
        "gnomAD_2": (
            "https://gnomad.broadinstitute.org/variant/"
            "{CHROM}-{POS}-{REF}-{ALT}?dataset=gnomad_r2_1"
        ),
        "autopvs1": "https://autopvs1.bgi.com/variant/hg19/{CHROM}-{POS}-{REF}-{ALT}",
        "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/?term={CHROM}-{POS}-{REF}-{ALT}",
    }

    # Find columns that contain URLs
    url_columns = []
    for col_idx, col_name in enumerate(df.columns):
        # Check if column contains URLs by sampling some values
        sample_values = df[col_name].dropna().astype(str).head(10)
        url_count = sum(
            1 for val in sample_values if val.startswith(("http://", "https://", "ftp://"))
        )

        if len(sample_values) > 0 and url_count / len(sample_values) > 0.5:
            url_columns.append((col_idx + 1, col_name))  # +1 because Excel is 1-indexed

    # Check if we have genomic coordinate columns for generating links
    has_genomic_coords = all(col in df.columns for col in ["CHROM", "POS", "REF", "ALT"])

    if has_genomic_coords:
        logger.info("Adding genomic link columns...")

        # Get the current max column to know where to add new columns
        start_col = ws.max_column + 1

        # Add headers for link columns
        for i, (link_name, template) in enumerate(link_templates.items()):
            col_letter = get_column_letter(start_col + i)
            ws[f"{col_letter}1"] = link_name

            # Add links for each data row
            for row_idx in range(2, ws.max_row + 1):
                # Get genomic coordinates from the current row
                chrom = (
                    ws[f"B{row_idx}"].value if "CHROM" in df.columns else None
                )  # Assuming CHROM is column B
                pos = (
                    ws[f"C{row_idx}"].value if "POS" in df.columns else None
                )  # Assuming POS is column C
                ref = (
                    ws[f"D{row_idx}"].value if "REF" in df.columns else None
                )  # Assuming REF is column D
                alt = (
                    ws[f"E{row_idx}"].value if "ALT" in df.columns else None
                )  # Assuming ALT is column E

                # Find actual column positions for genomic coordinates
                chrom_col = None
                pos_col = None
                ref_col = None
                alt_col = None

                for col_idx, col_name in enumerate(df.columns):
                    if col_name == "CHROM":
                        chrom_col = get_column_letter(col_idx + 1)
                        chrom = ws[f"{chrom_col}{row_idx}"].value
                    elif col_name == "POS":
                        pos_col = get_column_letter(col_idx + 1)
                        pos = ws[f"{pos_col}{row_idx}"].value
                    elif col_name == "REF":
                        ref_col = get_column_letter(col_idx + 1)
                        ref = ws[f"{ref_col}{row_idx}"].value
                    elif col_name == "ALT":
                        alt_col = get_column_letter(col_idx + 1)
                        alt = ws[f"{alt_col}{row_idx}"].value

                if all(x is not None for x in [chrom, pos, ref, alt]):
                    # Generate the URL
                    url = template.format(CHROM=chrom, POS=pos, REF=ref, ALT=alt)

                    # Set the cell with hyperlink
                    cell = ws[f"{col_letter}{row_idx}"]
                    cell.hyperlink = url
                    cell.value = "🔗"
                    cell.font = Font(color="0563C1", underline="single")

    # Add hyperlinks for existing URL columns
    for col_idx, col_name in url_columns:
        col_letter = get_column_letter(col_idx)
        logger.info(f"Processing existing links in column {col_name} ({col_letter})")

        # Iterate through data rows (starting from row 2, after header)
        for row_idx in range(2, ws.max_row + 1):
            cell = ws[f"{col_letter}{row_idx}"]
            if cell.value and str(cell.value).startswith(("http://", "https://", "ftp://")):
                # Set the hyperlink
                cell.hyperlink = str(cell.value)
                # Set display text to be more user-friendly
                cell.value = "🔗"
                # Apply blue color and underline (standard link formatting)
                cell.font = Font(color="0563C1", underline="single")

    # Update auto-filter to include new columns
    if has_genomic_coords:
        max_col_letter = get_column_letter(ws.max_column)
        ws.auto_filter.ref = f"A1:{max_col_letter}1"

    # Save the workbook
    wb.save(xlsx_file)
    genomic_links_count = len(link_templates) if has_genomic_coords else 0
    logger.info(
        f"Excel file finalized with {len(url_columns)} existing URL columns "
        f"and {genomic_links_count} genomic link columns"
    )


def create_report(output_dir, report_name, df, summary):
    """Create the HTML report using the template and data."""
    logger.info("Creating HTML report...")

    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Set up Jinja2 environment
    env = Environment(loader=FileSystemLoader(os.path.join(script_dir, "templates")))
    template = env.get_template("cohort_report_template.html")

    # Prepare column metadata for template
    column_metadata = prepare_column_metadata(df)

    # Convert DataFrame and summary to JSON strings for embedding in HTML
    variants_json = df.to_json(orient="records", date_format="iso")
    summary_json = json.dumps(summary)
    column_metadata_json = json.dumps(column_metadata)

    # Render the template with embedded JSON data
    html_content = template.render(
        report_name=report_name,
        generation_date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        variants_json=variants_json,
        summary_json=summary_json,
        column_metadata_json=column_metadata_json,
    )

    # Write the HTML to a file
    report_file = os.path.join(output_dir, "cohort_report.html")
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    logger.info(f"HTML report created: {report_file}")
    return report_file


def main():
    """Execute the cohort report creation."""
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

    # Create Excel report with all variants and clickable links
    excel_file = create_excel_report(df, args.output_dir)

    logger.info("Cohort report creation complete!")
    logger.info(f"HTML report available at: {report_file}")
    logger.info(f"Excel report available at: {excel_file}")


if __name__ == "__main__":
    main()

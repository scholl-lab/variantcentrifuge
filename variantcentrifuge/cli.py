# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

"""
CLI module for variantcentrifuge.

This module:
- Parses command line arguments.
- Sets up logging and configuration.
- Validates mandatory fields (reference, filters, and fields to extract) before
  proceeding with analysis, ensuring a clear error message is shown and the
  process exits early if these are missing.
- Validates input files (VCF, phenotype, and gene files) before processing.
  If required files or columns are missing, logs a clear error and exits
  gracefully, preventing partial execution.
- Invokes external tools like SnpSift, bcftools, and snpEff based on user inputs.
- Runs the analysis pipeline:
  1) Filters variants using SnpSift and bcftools.
  2) Extracts fields.
  3) Optionally replaces genotypes.
  4) Integrates phenotype data if provided.
  5) Runs variant-level analysis using analyze_variants (no gene burden) to
     produce a final TSV.
  6) If gene burden requested, runs analyze_variants again (with gene burden)
     to produce a separate TSV.
  7) Generates metadata, optionally converts to Excel, and appends additional
     sheets.

It ensures:
- The main output TSV (and XLSX if requested) always contains the variant-level
  results as the primary "Results" sheet.
- If gene burden is performed, its results appear in a separate "Gene Burden"
  sheet, leaving the main results intact.

Enhancement #13:
- Validate that mandatory fields (reference, filters, and fields) are provided
  either via CLI arguments or config file.
- If missing, log a clear error message and exit early.
- Prevent ambiguous behavior or silent failures when critical parameters are
  not set.

Enhancement #9:
- Validate input files (e.g., VCF, phenotype, gene files) before processing.
- Log clear errors and exit gracefully if required columns or data are missing.
- Avoid partial execution when prerequisites are not met.
"""

import argparse
import sys
import os
import hashlib
import datetime
import subprocess
import logging

from . import __version__
from .config import load_config
from .gene_bed import get_gene_bed
from .converter import convert_to_excel, append_tsv_as_sheet
from .phenotype_filter import filter_phenotypes
from .analyze_variants import analyze_variants
from .utils import check_external_tools, run_command, get_tool_version
from .replacer import replace_genotypes
from .phenotype import load_phenotypes, aggregate_phenotypes_for_samples

logger = logging.getLogger("variantcentrifuge")


def sanitize_metadata_field(value):
    """
    Sanitize a metadata field by removing tabs and newlines,
    replacing them with spaces to ensure TSV compatibility.

    Parameters
    ----------
    value : str
        The value to sanitize.

    Returns
    -------
    str
        Sanitized string value with no tabs or newlines.
    """
    return value.replace("\t", " ").replace("\n", " ").strip()


def normalize_genes(gene_name_str, gene_file_str):
    """
    Normalize genes from either a single gene name, a list of genes,
    or a file containing gene names.

    If 'all' is provided or no genes after filtering, returns "all".

    Parameters
    ----------
    gene_name_str : str or None
        The gene name(s) provided via CLI.
    gene_file_str : str or None
        Path to a file containing gene names.

    Returns
    -------
    str
        A normalized, space-separated string of gene names or "all".
    """
    if gene_file_str and gene_file_str.strip():
        if not os.path.exists(gene_file_str):
            logger.error(f"Gene file {gene_file_str} not found.")
            sys.exit(1)
        genes_from_file = []
        with open(gene_file_str, "r", encoding="utf-8") as gf:
            for line in gf:
                line = line.strip()
                if line:
                    genes_from_file.append(line)
        genes = genes_from_file
    else:
        if not gene_name_str:
            logger.error("No gene name provided and no gene file provided.")
            sys.exit(1)
        # Prevent confusion if a file was incorrectly given to -g
        if os.path.exists(gene_name_str):
            logger.error(
                f"It looks like you provided a file '{gene_name_str}' to "
                f"-g/--gene-name."
            )
            logger.error(
                "If you meant to provide a file of gene names, please use "
                "-G/--gene-file instead."
            )
            sys.exit(1)

        g_str = gene_name_str.replace(",", " ")
        genes = [
            g.strip()
            for g_str_part in g_str.split()
            for g in [g_str_part.strip()] if g
        ]

    if len(genes) == 1 and genes[0].lower() == "all":
        return "all"

    # Sort and deduplicate
    genes = sorted(set(genes))
    if not genes:
        return "all"
    return " ".join(genes)


def remove_vcf_extensions(filename):
    """
    Remove common VCF-related extensions from a filename.

    Parameters
    ----------
    filename : str
        The input filename, possibly ending in .vcf, .vcf.gz, or .gz

    Returns
    -------
    str
        The filename base without VCF-related extensions.
    """
    base = filename
    if base.endswith(".vcf.gz"):
        base = base[:-7]
    elif base.endswith(".vcf"):
        base = base[:-4]
    elif base.endswith(".gz"):
        base = base[:-3]
    return base


def compute_base_name(vcf_path, gene_name):
    """
    Compute a base name for output files based on the VCF filename and gene(s).

    If multiple genes are specified, create a hash to represent them.
    If 'all' is specified, append '.all'.
    Otherwise, append the gene name if it's not already in the VCF base name.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.
    gene_name : str
        The normalized gene name string.

    Returns
    -------
    str
        A base name for output files.
    """
    genes = gene_name.strip()
    vcf_base = os.path.basename(vcf_path)
    vcf_base = remove_vcf_extensions(vcf_base)

    if genes.lower() == "all":
        return f"{vcf_base}.all"
    split_genes = genes.split()
    if len(split_genes) > 1:
        gene_hash = hashlib.md5(genes.encode('utf-8')).hexdigest()[:8]
        return f"{vcf_base}.multiple-genes-{gene_hash}"
    else:
        if split_genes[0].lower() in vcf_base.lower():
            return vcf_base
        else:
            return f"{vcf_base}.{split_genes[0]}"


def parse_samples_from_vcf(vcf_file):
    """
    Parse samples from a VCF file by reading its header.

    Parameters
    ----------
    vcf_file : str
        Path to the VCF file.

    Returns
    -------
    list
        A list of sample names extracted from the VCF header.
    """
    output = run_command(["bcftools", "view", "-h", vcf_file], output_file=None)
    chrom_line = None
    for line in output.splitlines():
        if line.startswith("#CHROM"):
            chrom_line = line
            break
    if not chrom_line:
        return []
    fields = chrom_line.strip().split("\t")
    if len(fields) <= 9:
        return []
    samples = fields[9:]
    return samples


def load_terms_from_file(file_path):
    """
    Load terms (HPO terms, sample IDs, etc.) from a file, one per line.

    Parameters
    ----------
    file_path : str
        Path to a file containing one term per line.

    Returns
    -------
    list
        A list of terms loaded from the file.
    """
    terms = []
    if file_path:
        if not os.path.exists(file_path):
            logger.error(f"Required file not found: {file_path}")
            sys.exit(1)
        with open(file_path, "r", encoding="utf-8") as f:
            found_any = False
            for line in f:
                t = line.strip()
                if t:
                    found_any = True
                    terms.append(t)
            if not found_any:
                logger.error(f"File {file_path} is empty or invalid.")
                sys.exit(1)
    return terms


def validate_vcf_file(vcf_path):
    """
    Validate that the input VCF file exists and is readable.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file.

    Raises
    ------
    SystemExit
        If the VCF file is missing, not readable, or invalid.
    """
    if not vcf_path or not os.path.exists(vcf_path):
        logger.error(f"VCF file not found: {vcf_path}")
        sys.exit(1)
    if os.path.getsize(vcf_path) == 0:
        logger.error(f"VCF file {vcf_path} is empty.")
        sys.exit(1)


def validate_phenotype_file(phenotype_file, sample_col, value_col):
    """
    Validate phenotype file presence, non-empty, and required columns.

    Parameters
    ----------
    phenotype_file : str
        Path to phenotype file.
    sample_col : str
        Name of the sample column.
    value_col : str
        Name of the phenotype value column.

    Raises
    ------
    SystemExit
        If the phenotype file is missing, empty, or required columns are missing.
    """
    if not phenotype_file:
        # No phenotype file provided, no validation needed here.
        return

    if not os.path.exists(phenotype_file):
        logger.error(f"Phenotype file not found: {phenotype_file}")
        sys.exit(1)
    if os.path.getsize(phenotype_file) == 0:
        logger.error(f"Phenotype file {phenotype_file} is empty.")
        sys.exit(1)

    # Check columns by loading phenotypes
    # We'll load once here to validate, and if valid, load_phenotypes will be called again
    # or we can store the result if needed, but to keep code minimal, just validate now.
    with open(phenotype_file, "r", encoding="utf-8") as pf:
        header = pf.readline().strip()
        if not header:
            logger.error(f"Phenotype file {phenotype_file} has no header.")
            sys.exit(1)
        columns = header.split("\t") if "\t" in header else header.split(",")
        # Ensure required columns are present
        if sample_col not in columns:
            logger.error(
                f"Phenotype sample column '{sample_col}' not found in {phenotype_file}."
            )
            sys.exit(1)
        if value_col not in columns:
            logger.error(
                f"Phenotype value column '{value_col}' not found in {phenotype_file}."
            )
            sys.exit(1)
        # Check if at least one data line exists
        data_line = pf.readline()
        if not data_line.strip():
            logger.error(
                f"Phenotype file {phenotype_file} contains only a header and no data."
            )
            sys.exit(1)


def main():
    """
    Main entry point for variantcentrifuge CLI.

    Steps:
    - Parse arguments.
    - Configure logging and load config.
    - Validate that mandatory parameters (reference, filters, fields) are
      provided. If not, log error and exit.
    - Validate input files (VCF, phenotype, gene files, etc.) are present and
      have the necessary data and columns before proceeding.
      If validation fails, log error and exit.
    - Setup paths and run external tools (bcftools, SnpSift).
    - Filter and extract fields from VCF.
    - Optional genotype replacement and phenotype integration.
    - Run analyze_variants (no gene burden) to produce main variant-level
      results.
    - If gene burden requested, run analyze_variants again (with gene burden)
      to produce gene_burden_tsv.
    - Produce metadata, optionally convert to Excel and append sheets.
    - Cleanup intermediates if requested.
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

    # Phenotype-based grouping
    parser.add_argument("--case-phenotypes",
                        help="Comma-separated HPO terms defining case group")
    parser.add_argument("--control-phenotypes",
                        help="Comma-separated HPO terms defining control group")
    parser.add_argument("--case-phenotypes-file",
                        help="File with HPO terms for case group")
    parser.add_argument("--control-phenotypes-file",
                        help="File with HPO terms for control group")

    # Explicit case/control sample sets
    parser.add_argument("--case-samples",
                        help="Comma-separated sample IDs defining the case group")
    parser.add_argument("--control-samples",
                        help="Comma-separated sample IDs defining the control "
                             "group")
    parser.add_argument("--case-samples-file",
                        help="File with sample IDs for case group")
    parser.add_argument("--control-samples-file",
                        help="File with sample IDs for control group")

    # Gene burden mode and correction
    parser.add_argument("--gene-burden-mode", choices=["samples", "alleles"],
                        default="alleles",
                        help="Mode for gene burden calculation: 'samples' or "
                             "'alleles'")
    parser.add_argument("--correction-method", choices=["fdr", "bonferroni"],
                        default="fdr",
                        help="Multiple testing correction method for gene "
                             "burden test")

    args = parser.parse_args()

    # Map CLI log-level to logging
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

    # Validate gene inputs
    if args.gene_name and args.gene_file:
        logger.error(
            "You can provide either -g/--gene-name or -G/--gene-file, "
            "but not both."
        )
        sys.exit(1)
    if not args.gene_name and not args.gene_file:
        logger.error(
            "You must provide either a gene name using -g or a gene file "
            "using -G."
        )
        sys.exit(1)

    # Validate VCF file existence
    validate_vcf_file(args.vcf_file)

    check_external_tools()

    # Load configuration and reference
    cfg = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")

    # Validate mandatory fields, filters, and reference
    if not reference:
        logger.error(
            "A reference database must be specified either via --reference or "
            "in the configuration file. Please provide a valid reference."
        )
        sys.exit(1)

    if not filters:
        logger.error(
            "No filters provided. Filters must be specified either via "
            "--filters or in the configuration file."
        )
        sys.exit(1)

    if not fields:
        logger.error(
            "No fields to extract were provided. Fields must be specified "
            "either via --fields or in the configuration file."
        )
        sys.exit(1)

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
            args.phenotype_value_column
        )

    # Normalize genes
    gene_name = normalize_genes(
        args.gene_name if args.gene_name else "",
        args.gene_file
    )
    logger.debug(f"Normalized gene list: {gene_name}")

    # Determine final output file paths
    if args.output_file is not None:
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = os.path.join(
                args.output_dir, os.path.basename(args.output_file)
            )
    else:
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, gene_name)
        final_output = os.path.join(
            args.output_dir, f"{base_name}.final.tsv"
        )

    # Set run configuration in cfg
    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats
    cfg["gene_burden_mode"] = args.gene_burden_mode
    cfg["correction_method"] = args.correction_method

    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    # Determine stats output file if needed
    if not cfg["no_stats"] and not args.stats_output_file:
        cfg["stats_output_file"] = os.path.join(
            intermediate_dir, f"{base_name}.statistics.tsv"
        )
    else:
        cfg["stats_output_file"] = args.stats_output_file

    # Load phenotypes if provided
    phenotypes = {}
    use_phenotypes = False
    if (args.phenotype_file and args.phenotype_sample_column and
            args.phenotype_value_column):
        phenotypes = load_phenotypes(
            args.phenotype_file,
            args.phenotype_sample_column,
            args.phenotype_value_column
        )
        if not phenotypes:
            logger.error(
                f"No phenotype data loaded from {args.phenotype_file}. "
                "Check file formatting and column names."
            )
            sys.exit(1)
        use_phenotypes = True

    # Load phenotype terms and validate their files
    case_hpo_terms = []
    control_hpo_terms = []
    if args.case_phenotypes:
        case_hpo_terms = [
            t.strip() for t in args.case_phenotypes.split(",") if t.strip()
        ]
    case_hpo_terms += load_terms_from_file(args.case_phenotypes_file)

    if args.control_phenotypes:
        control_hpo_terms = [
            t.strip() for t in args.control_phenotypes.split(",") if t.strip()
        ]
    control_hpo_terms += load_terms_from_file(args.control_phenotypes_file)

    cfg["case_phenotypes"] = case_hpo_terms
    cfg["control_phenotypes"] = control_hpo_terms

    # Load explicit sample sets if provided
    case_samples = []
    control_samples = []
    if args.case_samples:
        case_samples = [
            s.strip() for s in args.case_samples.split(",") if s.strip()
        ]
    case_samples += load_terms_from_file(args.case_samples_file)

    if args.control_samples:
        control_samples = [
            s.strip() for s in args.control_samples.split(",") if s.strip()
        ]
    control_samples += load_terms_from_file(args.control_samples_file)

    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    # Extract gene regions
    logger.debug("Starting gene BED extraction.")
    bed_file = get_gene_bed(
        reference,
        gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True),
        output_dir=args.output_dir
    )
    logger.debug(f"Gene BED created at: {bed_file}")

    # Validate gene BED file
    if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
        logger.error(
            "Gene BED file could not be created or is empty. This may happen "
            "if no genes match your criteria or your reference is invalid."
        )
        sys.exit(1)

    # Check for missing genes
    if gene_name.lower() != "all":
        requested_genes = gene_name.split()
        found_genes = set()
        with open(bed_file, "r", encoding="utf-8") as bf:
            for line in bf:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    g_name = parts[3].split(";")[0].strip()
                    found_genes.add(g_name)
        missing = [g for g in requested_genes if g not in found_genes]
        if missing:
            logger.warning(
                "The following gene(s) were not found in the reference: "
                f"{', '.join(missing)}"
            )
            if cfg.get("debug_level", "INFO") == "ERROR":
                sys.exit(1)

    # Define intermediate and final paths
    variants_file = os.path.join(
        intermediate_dir, f"{base_name}.variants.vcf.gz"
    )
    filtered_file = os.path.join(
        intermediate_dir, f"{base_name}.filtered.vcf"
    )
    extracted_tsv = os.path.join(
        intermediate_dir, f"{base_name}.extracted.tsv"
    )
    genotype_replaced_tsv = os.path.join(
        intermediate_dir, f"{base_name}.genotype_replaced.tsv"
    )
    phenotype_added_tsv = os.path.join(
        intermediate_dir, f"{base_name}.phenotypes_added.tsv"
    )
    gene_burden_tsv = os.path.join(
        args.output_dir, f"{base_name}.gene_burden.tsv"
    )

    # Extract and filter variants
    logger.debug("Extracting variants and compressing as gzipped VCF...")
    run_command(["bcftools", "view", args.vcf_file, "-R", bed_file,
                 "-Oz", "-o", variants_file])
    logger.debug(f"Gzipped VCF saved to {variants_file}, indexing now...")
    run_command(["bcftools", "index", variants_file])
    logger.debug(f"Index created for {variants_file}")

    run_command(["SnpSift", "filter", filters, variants_file],
                output_file=filtered_file)
    logger.debug("Filter applied. Extracting fields...")

    field_list = fields.strip().split()
    run_command(["SnpSift", "extractFields", "-s", ",", "-e", "NA",
                 filtered_file] + field_list, output_file=extracted_tsv)
    if not os.path.exists(extracted_tsv) or os.path.getsize(extracted_tsv) == 0:
        logger.error(
            "Extraction of fields failed or produced empty results. Check "
            "your filters and fields configuration."
        )
        sys.exit(1)
    with open(extracted_tsv, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if lines:
        header = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")
        lines[0] = header
    with open(extracted_tsv, "w", encoding="utf-8") as f:
        f.writelines(lines)
    logger.debug(f"Fields extracted to {extracted_tsv}")

    # Check if GT column is present
    with open(extracted_tsv, "r", encoding="utf-8") as f:
        header_line = f.readline().strip()
    columns = header_line.split("\t")
    if "GT" not in columns:
        logger.warning(
            "GT column not found in extracted fields. No genotype data could "
            "be replaced. This may occur if your input VCF does not contain "
            "genotype data or if the 'GT' field was not included in the "
            "extraction fields. Ensure the input VCF contains genotypes and "
            "the 'GT' field is included in the extraction configuration."
        )
        if not args.no_replacement:
            logger.warning(
                "Genotype replacement was requested but cannot be performed "
                "due to missing 'GT' column. Your output will remain "
                "unchanged. Consider updating your input VCF or extraction "
                "parameters."
            )

    # Optional genotype replacement
    if not args.no_replacement and "GT" in columns:
        samples = parse_samples_from_vcf(variants_file)
        cfg["sample_list"] = ",".join(samples) if samples else ""
        logger.debug(
            f"Extracted {len(samples)} samples from VCF header for genotype "
            "replacement."
        )
        with open(extracted_tsv, "r", encoding="utf-8") as inp, \
                open(genotype_replaced_tsv, "w", encoding="utf-8") as out:
            for line in replace_genotypes(inp, cfg):
                out.write(line + "\n")
        replaced_tsv = genotype_replaced_tsv
        logger.debug(
            f"Genotype replacement done, results in {genotype_replaced_tsv}"
        )
    else:
        replaced_tsv = extracted_tsv
        logger.debug("Genotype replacement skipped or not applicable.")

    # Integrate phenotypes if provided
    if use_phenotypes:
        logger.debug("Adding phenotypes to the final table...")
        import re
        pattern = re.compile(r"^([^()]+)(?:\([^)]+\))?$")

        with open(replaced_tsv, "r", encoding="utf-8") as inp, \
                open(phenotype_added_tsv, "w", encoding="utf-8") as out:
            header = next(inp).rstrip("\n")
            header_fields = header.split("\t")
            header_fields.append("phenotypes")
            out.write("\t".join(header_fields) + "\n")

            gt_idx = header_fields.index("GT") if "GT" in header_fields else None

            wrote_data = False
            for line in inp:
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line + "\n")
                    continue
                fields_line = line.split("\t")

                if gt_idx is not None and gt_idx < len(fields_line):
                    gt_value = fields_line[gt_idx]
                    samples_in_line = []
                    if gt_value.strip():
                        sample_entries = gt_value.split(";")
                        for entry in sample_entries:
                            entry = entry.strip()
                            if entry:
                                m = pattern.match(entry)
                                if m:
                                    s = m.group(1).strip()
                                    if s:
                                        samples_in_line.append(s)
                    pheno_str = ""
                    if samples_in_line:
                        pheno_str = aggregate_phenotypes_for_samples(
                            samples_in_line, phenotypes
                        )
                    fields_line.append(pheno_str)
                else:
                    fields_line.append("")
                out.write("\t".join(fields_line) + "\n")
                wrote_data = True

            if not wrote_data:
                logger.error(
                    "Phenotype integration did not produce any output. Check "
                    "your phenotype data and parameters."
                )
                sys.exit(1)
        final_tsv = phenotype_added_tsv
    else:
        final_tsv = replaced_tsv

    # If Excel requested, set final_excel_file in cfg
    if args.xlsx:
        excel_file = os.path.splitext(final_output)[0] + ".xlsx"
        cfg["final_excel_file"] = excel_file

    # Produce the main variant-level results (no gene burden)
    temp_cfg = cfg.copy()
    temp_cfg["perform_gene_burden"] = False
    buffer = []
    with open(final_tsv, "r", encoding="utf-8") as inp:
        line_count = 0
        for line in analyze_variants(inp, temp_cfg):
            buffer.append(line)
            line_count += 1
        if line_count == 0:
            logger.error(
                "No variant-level results produced. Check your filters, "
                "fields, and input data."
            )
            sys.exit(1)

    if final_to_stdout:
        final_file = None
        if buffer:
            sys.stdout.write("\n".join(buffer) + "\n")
        final_out_path = None
    else:
        final_file = final_output
        with open(final_file, "w", encoding="utf-8") as out:
            for line in buffer:
                out.write(line + "\n")
        final_out_path = final_file

    # If gene burden is requested, produce separate gene_burden_tsv
    if cfg.get("perform_gene_burden", False):
        gene_burden_cfg = cfg.copy()
        gene_burden_cfg["perform_gene_burden"] = True
        line_count = 0
        with open(final_tsv, "r", encoding="utf-8") as inp, \
                open(gene_burden_tsv, "w", encoding="utf-8") as out:
            for line in analyze_variants(inp, gene_burden_cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count == 0:
            logger.warning(
                "Gene burden requested but no lines produced. "
                "No gene_burden_tsv created."
            )
        else:
            logger.debug(
                f"Gene burden analysis complete, results in {gene_burden_tsv}"
            )

    # Produce metadata
    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(
        f"Run ended at {end_time.isoformat()}, duration: {duration} seconds"
    )

    metadata_file = os.path.join(args.output_dir, f"{base_name}.metadata.tsv")
    with open(metadata_file, "w", encoding="utf-8") as mf:
        mf.write("Parameter\tValue\n")

        def meta_write(param, val):
            p = sanitize_metadata_field(param)
            v = sanitize_metadata_field(val)
            mf.write(f"{p}\t{v}\n")

        meta_write("Tool", "variantcentrifuge")
        meta_write("Version", __version__)
        meta_write("Run_start_time", start_time.isoformat())
        meta_write("Run_end_time", end_time.isoformat())
        meta_write("Run_duration_seconds", str(duration))
        meta_write("Date", datetime.datetime.now().isoformat())
        meta_write(
            "Command_line",
            " ".join([sanitize_metadata_field(x) for x in sys.argv])
        )

        for k, v in cfg.items():
            meta_write(f"config.{k}", str(v))

        # Unified version retrieval
        snpeff_ver = get_tool_version("snpEff")
        bcftools_ver = get_tool_version("bcftools")
        snpsift_ver = get_tool_version("SnpSift")
        bedtools_ver = get_tool_version("bedtools")

        meta_write("tool.snpeff_version", snpeff_ver)
        meta_write("tool.bcftools_version", bcftools_ver)
        meta_write("tool.snpsift_version", snpsift_ver)
        meta_write("tool.bedtools_version", bedtools_ver)

    # Convert to Excel if requested
    if args.xlsx and final_out_path and not final_to_stdout:
        if not os.path.exists(final_out_path) or \
                os.path.getsize(final_out_path) == 0:
            logger.warning(
                "Final output file is empty. Cannot convert to Excel."
            )
        else:
            logger.debug(
                "Converting final output (variant-level results) to Excel..."
            )
            xlsx_file = convert_to_excel(final_out_path, cfg)
            logger.debug("Excel conversion complete. Appending Metadata sheet...")
            append_tsv_as_sheet(xlsx_file, metadata_file, sheet_name="Metadata")

            if (not cfg["no_stats"] and cfg.get("stats_output_file") and
                    os.path.exists(cfg["stats_output_file"])):
                if os.path.getsize(cfg["stats_output_file"]) == 0:
                    logger.warning(
                        "Stats file is empty, skipping Statistics sheet."
                    )
                else:
                    logger.debug("Appending Statistics sheet to Excel...")
                    append_tsv_as_sheet(
                        xlsx_file, cfg["stats_output_file"],
                        sheet_name="Statistics"
                    )

            # If gene burden performed, append as separate sheet
            if (cfg.get("perform_gene_burden", False) and
                    os.path.exists(gene_burden_tsv) and
                    os.path.getsize(gene_burden_tsv) > 0):
                logger.debug(
                    "Appending Gene Burden results as a separate sheet..."
                )
                append_tsv_as_sheet(xlsx_file, gene_burden_tsv,
                                    sheet_name="Gene Burden")

    # Remove intermediates if requested
    if not args.keep_intermediates:
        logger.debug("Removing intermediate files...")
        intermediates = [variants_file, filtered_file, extracted_tsv]
        if not args.no_replacement and "GT" in columns:
            intermediates.append(genotype_replaced_tsv)
        if use_phenotypes:
            intermediates.append(phenotype_added_tsv)
        if final_file in intermediates:
            intermediates.remove(final_file)
        for f in intermediates:
            if os.path.exists(f):
                os.remove(f)
        logger.debug("Intermediate files removed.")

    logger.info(
        f"Processing completed. Output saved to "
        f"{'stdout' if final_to_stdout else final_output}"
    )


if __name__ == "__main__":
    main()

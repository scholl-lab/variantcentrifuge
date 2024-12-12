# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

"""
Command-line interface (CLI) module.

This script filters and processes VCF files using bcftools, snpEff, SnpSift, and optional steps.
Changes:
- --output-dir now defaults to "output".
- Intermediate files are placed in a "intermediate" subdirectory of --output-dir.
- Final outputs (final.tsv, final.xlsx, gene_burden.tsv if requested) remain in --output-dir.
- If --keep-intermediates is false, intermediate files are removed, but final outputs remain.

Other logic remains the same:
- --output-file optional, can be 'stdout', '-', or a filename. If not provided, compute base name.
- If gene burden is performed, produce a separate gene_burden file in output-dir.
- Convert to Excel if requested.

Requires other modules as previously discussed.
"""

import argparse
import sys
import os
import hashlib
from . import __version__
from .config import load_config
from .gene_bed import get_gene_bed
from .converter import convert_to_excel
from .replacer import replace_genotypes
from .phenotype_filter import filter_phenotypes
from .analyze_variants import analyze_variants
from .utils import log_message, check_external_tools, set_log_level, run_command

def compute_base_name(vcf_path, gene_name):
    vcf_basename = os.path.splitext(os.path.basename(vcf_path))[0]
    genes = gene_name.strip()
    if " " in genes or "," in genes:
        # multiple genes
        gene_hash = hashlib.md5(genes.encode('utf-8')).hexdigest()[:8]
        return f"{vcf_basename}.{gene_hash}"
    else:
        # single gene
        if genes.lower() in vcf_basename.lower():
            return vcf_basename
        else:
            return f"{vcf_basename}.{genes}"

def main():
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
    parser.add_argument("-c", "--config", help="Path to configuration file", default=None)
    parser.add_argument("-g", "--gene-name", help="Gene or list of genes of interest", required=True)
    parser.add_argument("-v", "--vcf-file", help="Input VCF file path", required=True)
    parser.add_argument("-r", "--reference", help="Reference database for snpEff", default=None)
    parser.add_argument("-f", "--filters", help="Filters to apply in SnpSift filter")
    parser.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    parser.add_argument("-s", "--samples-file", help="Samples file for genotype replacement")
    parser.add_argument("-o", "--output-file", nargs='?', const='stdout',
                        help="Final output file name or 'stdout'/'-' for stdout (if provided without argument, defaults to stdout)")
    parser.add_argument("--xlsx", action="store_true", help="Convert final result to Excel")
    parser.add_argument("--no-replacement", action="store_true", help="Skip genotype replacement step")
    parser.add_argument("--stats-output-file", help="File to write analysis statistics")
    parser.add_argument("--perform-gene-burden", action="store_true", help="Perform gene burden analysis")
    parser.add_argument("--output-dir", help="Directory to store intermediate and final output files", default="output")
    parser.add_argument("--keep-intermediates", action="store_true", default=True,
                        help="Keep intermediate files (default: true). If false, intermediate files are deleted at the end.")

    args = parser.parse_args()
    set_log_level(args.log_level)
    log_message("DEBUG", f"CLI arguments: {args}")

    check_external_tools()

    cfg = load_config(args.config)
    log_message("DEBUG", f"Configuration loaded: {cfg}")
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")

    if args.output_file is not None:
        # --output-file provided
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, args.gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = os.path.splitext(os.path.basename(args.output_file))[0]
            final_output = os.path.join(args.output_dir, os.path.basename(args.output_file))
    else:
        # No --output-file given
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, args.gene_name)
        final_output = os.path.join(args.output_dir, f"{base_name}.final.tsv")

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["stats_output_file"] = args.stats_output_file

    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    # BED file
    log_message("DEBUG", "Starting gene BED extraction.")
    bed_file = get_gene_bed(
        reference,
        args.gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True)
    )
    log_message("DEBUG", f"Gene BED created at: {bed_file}")

    # Intermediate files in intermediate_dir
    variants_file = os.path.join(intermediate_dir, f"{base_name}.variants.vcf")
    filtered_file = os.path.join(intermediate_dir, f"{base_name}.filtered.vcf")
    extracted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv")
    genotype_replaced_tsv = os.path.join(intermediate_dir, f"{base_name}.genotype_replaced.tsv")
    phenotype_filtered_tsv = os.path.join(intermediate_dir, f"{base_name}.phenotype_filtered.tsv")
    gene_burden_tsv = os.path.join(args.output_dir, f"{base_name}.gene_burden.tsv")  # gene burden in output-dir (final result)

    # Extract variants
    run_command(["bcftools", "view", args.vcf_file, "-R", bed_file], output_file=variants_file)
    log_message("DEBUG", f"Variants extracted to {variants_file}. Applying filters...")

    # Apply SnpSift filter
    run_command(["SnpSift", "filter", filters, variants_file], output_file=filtered_file)
    log_message("DEBUG", f"Filter applied. Extracting fields...")

    # Extract fields
    field_list = fields.strip().split()
    run_command(["SnpSift", "extractFields", "-s", ",", "-e", "NA", filtered_file] + field_list, output_file=extracted_tsv)
    # Clean header
    with open(extracted_tsv, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if lines:
        header = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")
        lines[0] = header
    with open(extracted_tsv, "w", encoding="utf-8") as f:
        f.writelines(lines)
    log_message("DEBUG", f"Fields extracted to {extracted_tsv}")

    # Genotype replacement
    if not args.no_replacement:
        log_message("DEBUG", "Replacing genotypes if sample info is available...")
        with open(extracted_tsv, "r", encoding="utf-8") as inp, open(genotype_replaced_tsv, "w", encoding="utf-8") as out:
            for line in replace_genotypes(inp, args.samples_file, cfg):
                out.write(line + "\n")
        replaced_tsv = genotype_replaced_tsv
        log_message("DEBUG", f"Genotype replacement done, results in {genotype_replaced_tsv}")
    else:
        replaced_tsv = extracted_tsv
        log_message("DEBUG", "Genotype replacement skipped.")

    # Phenotype filtering
    if cfg.get("use_phenotype_filtering", False):
        log_message("DEBUG", "Applying phenotype filtering...")
        with open(replaced_tsv, "r", encoding="utf-8") as inp, open(phenotype_filtered_tsv, "w", encoding="utf-8") as out:
            for line in filter_phenotypes(inp, cfg):
                out.write(line + "\n")
        replaced_tsv = phenotype_filtered_tsv
        log_message("DEBUG", "Phenotype filtering complete.")

    # Gene burden
    if cfg.get("perform_gene_burden", False):
        log_message("DEBUG", "Analyzing variants (gene burden) requested...")
        line_count = 0
        with open(replaced_tsv, "r", encoding="utf-8") as inp, open(gene_burden_tsv, "w", encoding="utf-8") as out:
            for line in analyze_variants(inp, cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count == 0:
            log_message("WARN", "No lines produced by analyze_variants. Falling back to replaced_tsv.")
            final_file = replaced_tsv
        else:
            final_file = replaced_tsv
        log_message("DEBUG", f"Variant analysis complete, gene burden results in {gene_burden_tsv}")
    else:
        log_message("DEBUG", "No gene burden analysis requested.")
        final_file = replaced_tsv

    # Write final output
    if final_to_stdout:
        log_message("DEBUG", "Writing final output to stdout.")
        with open(final_file, "r", encoding="utf-8") as inp:
            sys.stdout.write(inp.read())
        final_out_path = None
    else:
        log_message("DEBUG", f"Writing final output to {final_output}")
        with open(final_file, "r", encoding="utf-8") as inp, open(final_output, "w", encoding="utf-8") as out:
            for line in inp:
                out.write(line)
        final_out_path = final_output

    # Convert to Excel if requested and not stdout
    if args.xlsx and not final_to_stdout:
        log_message("DEBUG", "Converting final output to Excel...")
        convert_to_excel(final_out_path, cfg)
        log_message("DEBUG", "Excel conversion complete.")

    # Remove intermediates if keep_intermediates is false
    if not args.keep_intermediates:
        log_message("DEBUG", "Removing intermediate files...")
        intermediates = [variants_file, filtered_file, extracted_tsv]
        if not args.no_replacement:
            intermediates.append(genotype_replaced_tsv)
        if cfg.get("use_phenotype_filtering", False):
            intermediates.append(phenotype_filtered_tsv)
        # gene_burden_tsv is considered a final result (like stats), do not delete
        # final_file and final_output are final results, do not delete them
        # If final_file is in intermediates, remove it
        if final_file in intermediates:
            intermediates.remove(final_file)
        for f in intermediates:
            if os.path.exists(f):
                os.remove(f)
        # If intermediate_dir is empty, we might remove it or leave it. Let's leave it.
        log_message("DEBUG", "Intermediate files removed.")

    log_message("INFO", f"Processing completed. Output saved to {'stdout' if final_to_stdout else final_output}")

if __name__ == "__main__":
    main()

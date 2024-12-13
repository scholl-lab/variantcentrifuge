# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

"""
Command-line interface (CLI) module.

Changes:
- For SnpSift version retrieval, we now use subprocess directly with check=False
  to avoid errors if SnpSift returns a non-zero exit code.
- We then parse stdout from SnpSift's output and extract the version line if present.
"""

import argparse
import sys
import os
import hashlib
import datetime
import subprocess
from . import __version__
from .config import load_config
from .gene_bed import get_gene_bed
from .converter import convert_to_excel, append_tsv_as_sheet
from .phenotype_filter import filter_phenotypes
from .analyze_variants import analyze_variants
from .utils import log_message, check_external_tools, set_log_level, run_command
from .replacer import replace_genotypes
from .phenotype import load_phenotypes, aggregate_phenotypes_for_samples

def compute_base_name(vcf_path, gene_name):
    vcf_basename = os.path.splitext(os.path.basename(vcf_path))[0]
    genes = gene_name.strip()
    if " " in genes or "," in genes:
        # multiple genes
        gene_hash = hashlib.md5(genes.encode('utf-8')).hexdigest()[:8]
        return f"{vcf_basename}.{gene_hash}"
    else:
        if genes.lower() in vcf_basename.lower():
            return vcf_basename
        else:
            return f"{vcf_basename}.{genes}"

def parse_samples_from_vcf(vcf_file):
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

def sanitize_metadata_field(value):
    # Replace tabs and newlines with spaces
    return value.replace("\t", " ").replace("\n", " ").strip()

def safe_run_snpeff():
    try:
        o = run_command(["snpEff", "-version"], output_file=None)
        return o.strip().splitlines()[0] if o else "N/A"
    except:
        return "N/A"

def safe_run_bcftools():
    try:
        o = run_command(["bcftools", "--version"], output_file=None)
        return o.strip().splitlines()[0] if o else "N/A"
    except:
        return "N/A"

def safe_run_snpsift():
    # We'll directly run 'SnpSift annotate' using subprocess, ignoring errors
    # and parse stdout if available.
    try:
        result = subprocess.run(["SnpSift", "annotate"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                check=False)
        ver = "N/A"
        # Check stdout for a line starting with "SnpSift version"
        for line in result.stdout.splitlines():
            if line.startswith("SnpSift version"):
                ver = line
                break
        # If not found in stdout, maybe in stderr?
        if ver == "N/A":
            for line in result.stderr.splitlines():
                if line.startswith("SnpSift version"):
                    ver = line
                    break
        return ver
    except:
        return "N/A"

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
    parser.add_argument("-o", "--output-file", nargs='?', const='stdout',
                        help="Final output file name or 'stdout'/'-' for stdout (if provided without argument, defaults to stdout)")
    parser.add_argument("--xlsx", action="store_true", help="Convert final result to Excel")
    parser.add_argument("--no-replacement", action="store_true", help="Skip genotype replacement step")
    parser.add_argument("--stats-output-file", help="File to write analysis statistics")
    parser.add_argument("--perform-gene-burden", action="store_true", help="Perform gene burden analysis")
    parser.add_argument("--output-dir", help="Directory to store intermediate and final output files", default="output")
    parser.add_argument("--keep-intermediates", action="store_true", default=True,
                        help="Keep intermediate files (default: true). If false, intermediate files are deleted at the end.")

    # Phenotype arguments
    parser.add_argument("--phenotype-file", help="Path to phenotype file (.csv or .tsv)")
    parser.add_argument("--phenotype-sample-column", help="Name of the column containing sample IDs in phenotype file")
    parser.add_argument("--phenotype-value-column", help="Name of the column containing phenotype values in phenotype file")

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
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, args.gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = os.path.splitext(os.path.basename(args.output_file))[0]
            final_output = os.path.join(args.output_dir, os.path.basename(args.output_file))
    else:
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, args.gene_name)
        final_output = os.path.join(args.output_dir, f"{base_name}.final.tsv")

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["stats_output_file"] = args.stats_output_file

    # Phenotype logic
    phenotypes = {}
    use_phenotypes = False
    if args.phenotype_file and args.phenotype_sample_column and args.phenotype_value_column:
        phenotypes = load_phenotypes(args.phenotype_file, args.phenotype_sample_column, args.phenotype_value_column)
        use_phenotypes = True

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

    # Intermediate files
    variants_file = os.path.join(intermediate_dir, f"{base_name}.variants.vcf")
    filtered_file = os.path.join(intermediate_dir, f"{base_name}.filtered.vcf")
    extracted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv")
    genotype_replaced_tsv = os.path.join(intermediate_dir, f"{base_name}.genotype_replaced.tsv")
    phenotype_added_tsv = os.path.join(intermediate_dir, f"{base_name}.phenotypes_added.tsv")
    gene_burden_tsv = os.path.join(args.output_dir, f"{base_name}.gene_burden.tsv")

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
        samples = parse_samples_from_vcf(args.vcf_file)
        cfg["sample_list"] = ",".join(samples) if samples else ""
        log_message("DEBUG", f"Extracted {len(samples)} samples from VCF header for genotype replacement.")
        with open(extracted_tsv, "r", encoding="utf-8") as inp, open(genotype_replaced_tsv, "w", encoding="utf-8") as out:
            for line in replace_genotypes(inp, cfg):
                out.write(line + "\n")
        replaced_tsv = genotype_replaced_tsv
        log_message("DEBUG", f"Genotype replacement done, results in {genotype_replaced_tsv}")
    else:
        replaced_tsv = extracted_tsv
        log_message("DEBUG", "Genotype replacement skipped.")

    # Add phenotypes if available
    if use_phenotypes:
        log_message("DEBUG", "Adding phenotypes to the final table...")
        import re
        pattern = re.compile(r"^([^()]+)(?:\([^)]+\))?$")

        with open(replaced_tsv, "r", encoding="utf-8") as inp, open(phenotype_added_tsv, "w", encoding="utf-8") as out:
            header = next(inp).rstrip("\n")
            header_fields = header.split("\t")
            header_fields.append("phenotypes")
            out.write("\t".join(header_fields) + "\n")

            gt_idx = header_fields.index("GT") if "GT" in header_fields else None

            for line in inp:
                line = line.rstrip("\n")
                if not line.strip():
                    out.write(line+"\n")
                    continue
                fields = line.split("\t")

                if gt_idx is not None and gt_idx < len(fields):
                    gt_value = fields[gt_idx]
                    samples_in_line = []
                    if gt_value.strip():
                        sample_entries = gt_value.split(";")
                        for entry in sample_entries:
                            entry = entry.strip()
                            if not entry:
                                continue
                            m = pattern.match(entry)
                            if m:
                                s = m.group(1).strip()
                                if s:
                                    samples_in_line.append(s)
                    pheno_str = ""
                    if samples_in_line:
                        pheno_str = aggregate_phenotypes_for_samples(samples_in_line, phenotypes)
                    fields.append(pheno_str)
                else:
                    fields.append("")
                out.write("\t".join(fields) + "\n")
        final_tsv = phenotype_added_tsv
    else:
        final_tsv = replaced_tsv

    # Perform gene burden analysis if requested
    if cfg.get("perform_gene_burden", False):
        log_message("DEBUG", "Analyzing variants (gene burden) requested...")
        line_count = 0
        with open(final_tsv, "r", encoding="utf-8") as inp, open(gene_burden_tsv, "w", encoding="utf-8") as out:
            for line in analyze_variants(inp, cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count == 0:
            log_message("WARN", "No lines produced by analyze_variants. Falling back to final_tsv.")
            final_file = final_tsv
        else:
            final_file = final_tsv
        log_message("DEBUG", f"Variant analysis complete, gene burden results in {gene_burden_tsv}")
    else:
        log_message("DEBUG", "No gene burden analysis requested.")
        final_file = final_tsv

    if args.output_file is not None and args.output_file in ["stdout", "-"]:
        final_to_stdout = True

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

    # Create a metadata file as a TSV with two columns: Parameter and Value
    metadata_file = os.path.join(args.output_dir, f"{base_name}.metadata.tsv")
    with open(metadata_file, "w", encoding="utf-8") as mf:
        mf.write("Parameter\tValue\n")
        def meta_write(param, val):
            p = sanitize_metadata_field(param)
            v = sanitize_metadata_field(val)
            mf.write(f"{p}\t{v}\n")

        meta_write("Tool", "variantcentrifuge")
        meta_write("Version", __version__)
        meta_write("Date", datetime.datetime.now().isoformat())
        meta_write("Command_line", " ".join([sanitize_metadata_field(x) for x in sys.argv]))

        for k, v in cfg.items():
            meta_write(f"config.{k}", str(v))

        snpeff_ver = safe_run_snpeff()
        bcftools_ver = safe_run_bcftools()
        snpsift_ver = safe_run_snpsift()

        meta_write("tool.snpeff_version", snpeff_ver)
        meta_write("tool.bcftools_version", bcftools_ver)
        meta_write("tool.snpsift_version", snpsift_ver)

    # Convert to Excel if requested and not stdout and final_out_path is set
    if args.xlsx and not final_to_stdout and final_out_path:
        log_message("DEBUG", "Converting final output to Excel...")
        xlsx_file = convert_to_excel(final_out_path, cfg)
        log_message("DEBUG", "Excel conversion complete.")

        # Append Metadata sheet
        log_message("DEBUG", "Appending metadata as a sheet to Excel...")
        append_tsv_as_sheet(xlsx_file, metadata_file, sheet_name="Metadata")

    # Remove intermediates if keep_intermediates is false
    if not args.keep_intermediates:
        log_message("DEBUG", "Removing intermediate files...")
        intermediates = [variants_file, filtered_file, extracted_tsv]
        if not args.no_replacement:
            intermediates.append(genotype_replaced_tsv)
        if use_phenotypes:
            intermediates.append(phenotype_added_tsv)
        # gene_burden_tsv and metadata_file are final results, do not remove
        # final_file and final_output are final results
        if final_file in intermediates:
            intermediates.remove(final_file)
        for f in intermediates:
            if os.path.exists(f):
                os.remove(f)
        log_message("DEBUG", "Intermediate files removed.")

    log_message("INFO", f"Processing completed. Output saved to {'stdout' if final_to_stdout else final_output}")

if __name__ == "__main__":
    main()

# File: variantcentrifuge/cli.py
# Location: variantcentrifuge/variantcentrifuge/cli.py

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
from .utils import check_external_tools, run_command
from .replacer import replace_genotypes
from .phenotype import load_phenotypes, aggregate_phenotypes_for_samples

logger = logging.getLogger("variantcentrifuge")

def sanitize_metadata_field(value):
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
    try:
        result = subprocess.run(["SnpSift", "annotate"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                text=True,
                                check=False)
        ver = "N/A"
        for line in result.stdout.splitlines():
            if line.startswith("SnpSift version"):
                ver = line
                break
        if ver == "N/A":
            for line in result.stderr.splitlines():
                if line.startswith("SnpSift version"):
                    ver = line
                    break
        return ver
    except:
        return "N/A"

def normalize_genes(gene_name_str, gene_file_str):
    if gene_file_str and gene_file_str.strip():
        if not os.path.exists(gene_file_str):
            logger.error(f"Gene file {gene_file_str} not found.")
            sys.exit(1)
        genes_from_file = []
        with open(gene_file_str, "r", encoding="utf-8") as gf:
            for line in gf:
                line=line.strip()
                if line:
                    genes_from_file.append(line)
        genes = genes_from_file
    else:
        if not gene_name_str:
            logger.error("No gene name provided and no gene file provided.")
            sys.exit(1)
        # Check if user gave a file to -g
        if os.path.exists(gene_name_str):
            logger.error(f"It looks like you provided a file '{gene_name_str}' to -g/--gene-name.")
            logger.error("If you meant to provide a file of gene names, please use -G/--gene-file instead.")
            sys.exit(1)

        g_str = gene_name_str.replace(",", " ")
        genes = [g.strip() for g_str_part in g_str.split() for g in [g_str_part.strip()] if g]

    if len(genes) == 1 and genes[0].lower() == "all":
        return "all"
    # sort and deduplicate
    genes = sorted(set(genes))
    if not genes:
        return "all"
    return " ".join(genes)

def remove_vcf_extensions(filename):
    base = filename
    if base.endswith(".vcf.gz"):
        base = base[:-7]
    elif base.endswith(".vcf"):
        base = base[:-4]
    elif base.endswith(".gz"):
        base = base[:-3]
    return base

def compute_base_name(vcf_path, gene_name):
    genes = gene_name.strip()
    vcf_base = os.path.basename(vcf_path)
    vcf_base = remove_vcf_extensions(vcf_base)

    if genes.lower() == "all":
        return f"{vcf_base}.all"
    else:
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
    """Load HPO terms or sample IDs from a file, one per line."""
    terms = []
    if file_path and os.path.exists(file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                t = line.strip()
                if t:
                    terms.append(t)
    return terms

def main():
    # Setup basic logging. Later, level will be set based on CLI argument.
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)

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
    parser.add_argument("-c", "--config", help="Path to configuration file", default=None)
    parser.add_argument("-g", "--gene-name", help="Gene or list of genes of interest")
    parser.add_argument("-G", "--gene-file", help="File containing gene names, one per line")
    parser.add_argument("-v", "--vcf-file", help="Input VCF file path", required=True)
    parser.add_argument("-r", "--reference", help="Reference database for snpEff", default=None)
    parser.add_argument("-f", "--filters", help="Filters to apply in SnpSift filter")
    parser.add_argument("-e", "--fields", help="Fields to extract with SnpSift extractFields")
    parser.add_argument("-o", "--output-file", nargs='?', const='stdout',
                        help="Final output file name or 'stdout'/'-' for stdout")
    parser.add_argument("--xlsx", action="store_true", help="Convert final result to Excel")
    parser.add_argument("--no-replacement", action="store_true", help="Skip genotype replacement step")
    parser.add_argument("--stats-output-file", help="File to write analysis statistics")
    parser.add_argument("--perform-gene-burden", action="store_true", help="Perform gene burden analysis")
    parser.add_argument("--output-dir", help="Directory to store intermediate and final output files", default="output")
    parser.add_argument("--keep-intermediates", action="store_true", default=True,
                        help="Keep intermediate files.")
    parser.add_argument("--phenotype-file", help="Path to phenotype file (.csv or .tsv)")
    parser.add_argument("--phenotype-sample-column", help="Name of the column containing sample IDs in phenotype file")
    parser.add_argument("--phenotype-value-column", help="Name of the column containing phenotype values in phenotype file")
    parser.add_argument("--no-stats", action="store_true", default=False,
                        help="Skip the statistics computation step.")

    # New CLI options for phenotype-driven grouping
    parser.add_argument("--case-phenotypes", help="Comma-separated HPO terms defining case group")
    parser.add_argument("--control-phenotypes", help="Comma-separated HPO terms defining control group")
    parser.add_argument("--case-phenotypes-file", help="File with HPO terms for case group")
    parser.add_argument("--control-phenotypes-file", help="File with HPO terms for control group")

    # New CLI options for explicit case/control sample specification
    parser.add_argument("--case-samples", help="Comma-separated sample IDs defining the case group")
    parser.add_argument("--control-samples", help="Comma-separated sample IDs defining the control group")
    parser.add_argument("--case-samples-file", help="File with sample IDs for case group")
    parser.add_argument("--control-samples-file", help="File with sample IDs for control group")

    args = parser.parse_args()

    # Map CLI log-level to logging module
    log_level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARN": logging.WARNING,
        "ERROR": logging.ERROR
    }
    logging.getLogger("variantcentrifuge").setLevel(log_level_map[args.log_level])

    logger.debug(f"CLI arguments: {args}")

    if args.gene_name and args.gene_file:
        logger.error("You can provide either -g/--gene-name or -G/--gene-file, but not both.")
        sys.exit(1)
    if not args.gene_name and not args.gene_file:
        logger.error("You must provide either a gene name using -g or a gene file using -G.")
        sys.exit(1)

    check_external_tools()

    cfg = load_config(args.config)
    logger.debug(f"Configuration loaded: {cfg}")
    reference = args.reference or cfg.get("reference")
    filters = args.filters or cfg.get("filters")
    fields = args.fields or cfg.get("fields_to_extract")

    gene_name = normalize_genes(args.gene_name if args.gene_name else "", args.gene_file)
    logger.debug(f"Normalized gene list: {gene_name}")

    if args.output_file is not None:
        if args.output_file in ["stdout", "-"]:
            final_to_stdout = True
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = None
        else:
            final_to_stdout = False
            base_name = compute_base_name(args.vcf_file, gene_name)
            final_output = os.path.join(args.output_dir, os.path.basename(args.output_file))
    else:
        final_to_stdout = False
        base_name = compute_base_name(args.vcf_file, gene_name)
        final_output = os.path.join(args.output_dir, f"{base_name}.final.tsv")

    cfg["perform_gene_burden"] = args.perform_gene_burden
    cfg["no_stats"] = args.no_stats

    os.makedirs(args.output_dir, exist_ok=True)
    intermediate_dir = os.path.join(args.output_dir, "intermediate")
    os.makedirs(intermediate_dir, exist_ok=True)

    if not cfg["no_stats"] and not args.stats_output_file:
        cfg["stats_output_file"] = os.path.join(intermediate_dir, f"{base_name}.statistics.tsv")
    else:
        cfg["stats_output_file"] = args.stats_output_file

    phenotypes = {}
    use_phenotypes = False
    if args.phenotype_file and args.phenotype_sample_column and args.phenotype_value_column:
        phenotypes = load_phenotypes(args.phenotype_file, args.phenotype_sample_column, args.phenotype_value_column)
        use_phenotypes = True

    # Load phenotype terms for case and control
    case_hpo_terms = []
    control_hpo_terms = []
    if args.case_phenotypes:
        case_hpo_terms = [t.strip() for t in args.case_phenotypes.split(",") if t.strip()]
    case_hpo_terms += load_terms_from_file(args.case_phenotypes_file)

    if args.control_phenotypes:
        control_hpo_terms = [t.strip() for t in args.control_phenotypes.split(",") if t.strip()]
    control_hpo_terms += load_terms_from_file(args.control_phenotypes_file)

    cfg["case_phenotypes"] = case_hpo_terms
    cfg["control_phenotypes"] = control_hpo_terms

    # Load explicit case/control sample sets if provided
    case_samples = []
    control_samples = []
    if args.case_samples:
        case_samples = [s.strip() for s in args.case_samples.split(",") if s.strip()]
    case_samples += load_terms_from_file(args.case_samples_file)

    if args.control_samples:
        control_samples = [s.strip() for s in args.control_samples.split(",") if s.strip()]
    control_samples += load_terms_from_file(args.control_samples_file)

    cfg["case_samples"] = case_samples
    cfg["control_samples"] = control_samples

    logger.debug("Starting gene BED extraction.")
    bed_file = get_gene_bed(
        reference,
        gene_name,
        interval_expand=cfg.get("interval_expand", 0),
        add_chr=cfg.get("add_chr", True),
        output_dir=args.output_dir
    )
    logger.debug(f"Gene BED created at: {bed_file}")

    if gene_name.lower() != "all":
        requested_genes = gene_name.split()
        found_genes = set()
        with open(bed_file, "r", encoding="utf-8") as bf:
            for line in bf:
                parts=line.strip().split("\t")
                if len(parts)>=4:
                    g_name = parts[3].split(";")[0].strip()
                    found_genes.add(g_name)
        missing = [g for g in requested_genes if g not in found_genes]
        if missing:
            logger.warning(f"The following gene(s) were not found in the reference: {', '.join(missing)}")
            if cfg.get("debug_level","INFO") == "ERROR":
                sys.exit(1)

    variants_file = os.path.join(intermediate_dir, f"{base_name}.variants.vcf.gz")
    filtered_file = os.path.join(intermediate_dir, f"{base_name}.filtered.vcf")
    extracted_tsv = os.path.join(intermediate_dir, f"{base_name}.extracted.tsv")
    genotype_replaced_tsv = os.path.join(intermediate_dir, f"{base_name}.genotype_replaced.tsv")
    phenotype_added_tsv = os.path.join(intermediate_dir, f"{base_name}.phenotypes_added.tsv")
    gene_burden_tsv = os.path.join(args.output_dir, f"{base_name}.gene_burden.tsv")

    logger.debug("Extracting variants and compressing as gzipped VCF with bcftools view...")
    run_command(["bcftools", "view", args.vcf_file, "-R", bed_file, "-Oz", "-o", variants_file])
    logger.debug(f"Gzipped VCF saved to {variants_file}, indexing now...")
    run_command(["bcftools", "index", variants_file])
    logger.debug(f"Index created for {variants_file}")

    run_command(["SnpSift", "filter", filters, variants_file], output_file=filtered_file)
    logger.debug("Filter applied. Extracting fields...")

    field_list = fields.strip().split()
    run_command(["SnpSift", "extractFields", "-s", ",", "-e", "NA", filtered_file] + field_list, output_file=extracted_tsv)
    with open(extracted_tsv, "r", encoding="utf-8") as f:
        lines = f.readlines()
    if lines:
        header = lines[0].replace("ANN[0].", "").replace("GEN[*].", "")
        lines[0] = header
    with open(extracted_tsv, "w", encoding="utf-8") as f:
        f.writelines(lines)
    logger.debug(f"Fields extracted to {extracted_tsv}")

    if not args.no_replacement:
        samples = parse_samples_from_vcf(variants_file)
        cfg["sample_list"] = ",".join(samples) if samples else ""
        logger.debug(f"Extracted {len(samples)} samples from VCF header for genotype replacement.")
        with open(extracted_tsv, "r", encoding="utf-8") as inp, open(genotype_replaced_tsv, "w", encoding="utf-8") as out:
            for line in replace_genotypes(inp, cfg):
                out.write(line + "\n")
        replaced_tsv = genotype_replaced_tsv
        logger.debug(f"Genotype replacement done, results in {genotype_replaced_tsv}")
    else:
        replaced_tsv = extracted_tsv
        logger.debug("Genotype replacement skipped.")

    if use_phenotypes:
        logger.debug("Adding phenotypes to the final table...")
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

    # Pass configuration down - analyze_variants will handle case/control assignment
    if cfg.get("perform_gene_burden", False):
        logger.debug("Analyzing variants (gene burden) requested...")
        line_count = 0
        with open(final_tsv, "r", encoding="utf-8") as inp, open(gene_burden_tsv, "w", encoding="utf-8") as out:
            for line in analyze_variants(inp, cfg):
                out.write(line + "\n")
                line_count += 1
        if line_count == 0:
            logger.warning("No lines produced by analyze_variants. Falling back to final_tsv.")
            final_file = final_tsv
        else:
            final_file = final_tsv
        logger.debug(f"Variant analysis complete, gene burden results in {gene_burden_tsv}")
    else:
        logger.debug("No gene burden analysis requested.")
        buffer = []
        with open(final_tsv, "r", encoding="utf-8") as inp:
            for line in analyze_variants(inp, cfg):
                buffer.append(line)
        if final_to_stdout:
            final_file = None
        else:
            final_file = final_output
            with open(final_file, "w", encoding="utf-8") as out:
                for line in buffer:
                    out.write(line + "\n")

    if args.output_file is not None and args.output_file in ["stdout", "-"]:
        final_to_stdout = True

    if final_to_stdout:
        logger.debug("Writing final output to stdout.")
        if final_file and os.path.exists(final_file):
            with open(final_file, "r", encoding="utf-8") as inp:
                sys.stdout.write(inp.read())
        final_out_path = None
    else:
        final_out_path = final_file

    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    logger.info(f"Run ended at {end_time.isoformat()}, duration: {duration} seconds")

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
        meta_write("Command_line", " ".join([sanitize_metadata_field(x) for x in sys.argv]))

        for k, v in cfg.items():
            meta_write(f"config.{k}", str(v))

        snpeff_ver = safe_run_snpeff()
        bcftools_ver = safe_run_bcftools()
        snpsift_ver = safe_run_snpsift()

        meta_write("tool.snpeff_version", snpeff_ver)
        meta_write("tool.bcftools_version", bcftools_ver)
        meta_write("tool.snpsift_version", snpsift_ver)

    if args.xlsx and final_out_path and not final_to_stdout:
        if not os.path.exists(final_out_path) or os.path.getsize(final_out_path) == 0:
            logger.warning("Final output file is empty. Cannot convert to Excel.")
        else:
            logger.debug("Converting final output to Excel...")
            xlsx_file = convert_to_excel(final_out_path, cfg)
            logger.debug("Excel conversion complete.")

            logger.debug("Appending metadata as a sheet to Excel...")
            append_tsv_as_sheet(xlsx_file, metadata_file, sheet_name="Metadata")

            if not cfg["no_stats"] and cfg.get("stats_output_file") and os.path.exists(cfg["stats_output_file"]):
                if os.path.getsize(cfg["stats_output_file"]) == 0:
                    logger.warning("Stats file is empty, skipping Statistics sheet.")
                else:
                    logger.debug("Appending statistics as a sheet to Excel...")
                    append_tsv_as_sheet(xlsx_file, cfg["stats_output_file"], sheet_name="Statistics")

    if not args.keep_intermediates:
        logger.debug("Removing intermediate files...")
        intermediates = [variants_file, filtered_file, extracted_tsv]
        if not args.no_replacement:
            intermediates.append(genotype_replaced_tsv)
        if use_phenotypes:
            intermediates.append(phenotype_added_tsv)
        if final_file in intermediates:
            intermediates.remove(final_file)
        for f in intermediates:
            if os.path.exists(f):
                os.remove(f)
        logger.debug("Intermediate files removed.")

    logger.info(f"Processing completed. Output saved to {'stdout' if final_to_stdout else final_output}")

if __name__ == "__main__":
    main()

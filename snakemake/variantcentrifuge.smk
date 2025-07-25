# snakemake/variantcentrifuge.smk
import os
import glob
import functools

# --- Configuration ---
configfile: "config_vc.yaml" # Assumes config_vc.yaml is in the same directory or specified via --configfile

# --- Global Parameters from Config ---
VCF_LIST_FILE = config["vcf_list_file"]
BASE_OUTPUT_FOLDER = config["base_output_folder"]
LOG_SUBFOLDER = config.get("log_subfolder", "logs_vc") # Default if not in config

VC_CONFIG_FILE = config["variantcentrifuge_config_file"]
VC_GENES = config["genes_of_interest"]
VC_LOG_LEVEL = config["log_level"]
VC_THREADS = config["threads_per_job"]
VC_XLSX = config.get("generate_xlsx", False)
VC_HTML_REPORT = config.get("generate_html_report", False)
VC_ADD_CHR = config.get("add_chr_prefix", False)
VC_PRESETS = config.get("presets", [])
VC_APPEND_GT_FIELDS = config.get("append_genotype_fields", [])
VC_GENE_LIST_FILES = config.get("gene_list_files", [])
VC_ENABLE_IGV = config.get("enable_igv", False)
VC_IGV_REF = config.get("igv_reference", "")
VC_BAM_MAPPING = config.get("bam_mapping_file", "")
VC_IGV_FASTA = config.get("igv_fasta", "")
VC_IGV_IDEOGRAM = config.get("igv_ideogram", "")
# MODIFIED: Start of IGV flanking feature
VC_IGV_FLANKING = config.get("igv_flanking", 50)
# MODIFIED: End of IGV flanking feature

CONDA_ENV_VC = config.get("conda_environment_variantcentrifuge", None)

# --- HPC Resource Defaults ---
DEFAULT_MEM_MB = config.get("default_job_mem_mb", 16000)
DEFAULT_TIME = config.get("default_job_time", "24:00:00")
SCRATCH_DIR = os.environ.get('TMPDIR', '/tmp') # Use system TMPDIR or HPC $TMPDIR

# --- Helper Functions ---
def get_vcf_info():
    """Load VCF file paths from a list file and extract sample information.
    Returns a dictionary mapping sample IDs to their full VCF paths.
    """
    if not os.path.exists(VCF_LIST_FILE):
        raise FileNotFoundError(f"VCF list file not found: {VCF_LIST_FILE}")

    vcf_paths = {}
    with open(VCF_LIST_FILE, 'r') as f:
        for line in f:
            path = line.strip()
            if path and os.path.exists(path):
                # Extract sample name from the file name (without directory and extension)
                sample = os.path.basename(path)
                if sample.endswith(".vcf.gz"):
                    sample = sample[:-7]  # Remove .vcf.gz
                elif sample.endswith(".vcf"):
                    sample = sample[:-4]  # Remove .vcf
                vcf_paths[sample] = path

    if not vcf_paths:
        raise ValueError(f"No valid VCF paths found in {VCF_LIST_FILE}")

    return vcf_paths

# Map of sample IDs to their VCF paths
VCF_PATHS = get_vcf_info()
SAMPLES = list(VCF_PATHS.keys())

def get_sample_output_dir(wildcards):
    """Define a unique output directory for each sample within the base output folder."""
    return os.path.join(BASE_OUTPUT_FOLDER, wildcards.sample)

def get_final_tsv_output(wildcards):
    """Determine the final TSV output name based on XLSX flag."""
    sample_out_dir = get_sample_output_dir(wildcards)
    # VariantCentrifuge typically names its output based on input VCF and gene.
    # For simplicity in Snakemake, we'll define a consistent output name.
    # The actual name chosen by variantcentrifuge inside its output_dir might differ slightly,
    # but we need *a* target for Snakemake. Let's assume variantcentrifuge --output-file controls this.
    return os.path.join(sample_out_dir, f"{wildcards.sample}.vc_analysis.tsv")

def get_final_xlsx_output(wildcards):
    """Define the final XLSX output name if requested."""
    sample_out_dir = get_sample_output_dir(wildcards)
    return os.path.join(sample_out_dir, f"{wildcards.sample}.vc_analysis.xlsx")

def get_html_report_output(wildcards):
    """Define the HTML report index file as a target if requested."""
    sample_out_dir = get_sample_output_dir(wildcards)
    return os.path.join(sample_out_dir, "report", "index.html")

# --- Target Rule ---
rule all:
    input:
        expand(os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "{sample}.vc_analysis.tsv"), sample=SAMPLES),
        expand(os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "{sample}.vc_analysis.xlsx"), sample=SAMPLES) if VC_XLSX else [],
        expand(os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "report", "index.html"), sample=SAMPLES) if VC_HTML_REPORT else []

# --- Main VariantCentrifuge Rule ---
rule run_variantcentrifuge:
    input:
        # Use the exact VCF path from our mapping
        vcf=lambda wildcards: VCF_PATHS[wildcards.sample]
    output:
        # Snakemake needs to track specific output files.
        # VariantCentrifuge's --output-dir will contain many files.
        # We track the main TSV, and conditionally the XLSX and HTML report.
        final_tsv = os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "{sample}.vc_analysis.tsv"),
        final_xlsx = os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "{sample}.vc_analysis.xlsx") if VC_XLSX else [], # Optional output
        html_report_index = os.path.join(BASE_OUTPUT_FOLDER, "{sample}", "report", "index.html") if VC_HTML_REPORT else []
    params:
        sample_output_dir = get_sample_output_dir,
        vc_config = VC_CONFIG_FILE,
        genes = VC_GENES,
        log_level = VC_LOG_LEVEL,
        threads = VC_THREADS,
        add_chr = "--add-chr" if VC_ADD_CHR else "",
        xlsx_flag = "--xlsx" if VC_XLSX else "",
        html_flag = "--html-report" if VC_HTML_REPORT else "",
        presets_flags = " ".join([f"--preset {p}" for p in VC_PRESETS]),
        append_gt_flags = ("--append-extra-sample-fields " + " ".join(VC_APPEND_GT_FIELDS)) if VC_APPEND_GT_FIELDS else "",
        gene_list_flags = " ".join([f"--annotate-gene-list {path}" for path in VC_GENE_LIST_FILES]) if VC_GENE_LIST_FILES else "",
        igv_flag = "--igv" if VC_ENABLE_IGV else "",
        igv_ref = f"--igv-reference {VC_IGV_REF}" if VC_ENABLE_IGV and VC_IGV_REF else "",
        bam_map = f"--bam-mapping-file {VC_BAM_MAPPING}" if VC_ENABLE_IGV and VC_BAM_MAPPING else "",
        igv_fasta = f"--igv-fasta {VC_IGV_FASTA}" if VC_ENABLE_IGV and VC_IGV_FASTA else "",
        igv_ideogram = f"--igv-ideogram {VC_IGV_IDEOGRAM}" if VC_ENABLE_IGV and VC_IGV_IDEOGRAM else "",
        # MODIFIED: Start of IGV flanking feature
        igv_flanking = f"--igv-flanking {VC_IGV_FLANKING}" if VC_ENABLE_IGV else ""
        # MODIFIED: End of IGV flanking feature
    log:
        os.path.join(BASE_OUTPUT_FOLDER, LOG_SUBFOLDER, "{sample}.variantcentrifuge.log")
    conda:
        CONDA_ENV_VC if CONDA_ENV_VC else None
    threads:
        VC_THREADS
    resources:
        mem_mb=DEFAULT_MEM_MB,
        time=DEFAULT_TIME,
        tmpdir=SCRATCH_DIR
    shell:
        """
        set -e -o pipefail # Fail on error and on undefined variables in pipes

        echo "Starting VariantCentrifuge for sample {wildcards.sample} at: $(date)" > {log}
        echo "Output directory: {params.sample_output_dir}" >> {log}

        # Create sample-specific output directory for variantcentrifuge
        # This is where variantcentrifuge will write its 'intermediate', 'report' subdirs etc.
        mkdir -p {params.sample_output_dir}

        variantcentrifuge \
            -c {params.vc_config} \
            -g {params.genes} \
            -v {input.vcf} \
            --output-file {output.final_tsv} \
            --output-dir {params.sample_output_dir} \
            --log-level {params.log_level} \
            --threads {threads} \
            {params.add_chr} \
            {params.xlsx_flag} \
            {params.html_flag} \
            {params.presets_flags} \
            {params.append_gt_flags} \
            {params.gene_list_flags} \
            {params.igv_flag} \
            {params.igv_ref} \
            {params.bam_map} \
            {params.igv_fasta} \
            {params.igv_ideogram} \
            {params.igv_flanking} \
            >> {log} 2>&1

        echo "Finished VariantCentrifuge for sample {wildcards.sample} at: $(date)" >> {log}
        """

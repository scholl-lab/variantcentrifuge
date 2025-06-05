# VariantCentrifuge Snakemake Workflow

This directory contains a Snakemake workflow designed to automate the execution of the `variantcentrifuge` tool across multiple VCF files.

## Overview

The workflow processes multiple VCF files (specified in a list file) using `variantcentrifuge`, applying a common set of parameters defined in a configuration file. It manages output directories, logging, and can be submitted to an HPC cluster using the provided SLURM script.

## Files

-   `variantcentrifuge.smk`: The main Snakemake workflow definition file.
-   `config_vc.yaml`: An example YAML configuration file. **You must copy and adapt this file for your specific paths and parameters.**
-   `run_variantcentrifuge.sh`: A Bash script for submitting the Snakemake workflow to a SLURM-managed HPC cluster.
-   `README_snakemake.md`: This file.

## Prerequisites

1.  **Snakemake:** Ensure Snakemake is installed (e.g., via conda or pip).
2.  **Conda:** If using conda environments for dependencies (recommended), ensure conda is installed.
3.  **VariantCentrifuge:** The `variantcentrifuge` command-line tool must be installed and accessible in the environment where Snakemake jobs will run, or a conda environment for it must be specified in `config_vc.yaml`.
4.  **External Tools for VariantCentrifuge:** All external tools required by `variantcentrifuge` (bcftools, SnpEff, SnpSift, bedtools, igv-reports) must be available in the execution environment or conda environment.
5.  **SLURM (Optional):** If running on an HPC, SLURM should be available, and you might need a Snakemake SLURM profile.

## Setup and Configuration

1.  **Copy and Edit Configuration:**
    Make a copy of `config_vc.yaml` (e.g., `my_vc_config.yaml`) and edit it to reflect your setup:
    *   `vcf_list_file`: Path to a text file listing the full paths to VCF files you want to process (one path per line).
    *   `base_output_folder`: Where all results will be written. A subdirectory will be created for each sample.
    *   `variantcentrifuge_config_file`: Path to the main JSON configuration file used by `variantcentrifuge` itself.
    *   Update other `variantcentrifuge` parameters (`genes_of_interest`, `threads_per_job`, flags, presets, IGV settings, etc.) as needed.
    *   `conda_environment_variantcentrifuge`: Specify the name of the conda environment that contains `variantcentrifuge` and its dependencies if you are using `--use-conda`.

2.  **Input VCFs List File:**
    Create a text file (e.g., `vcf_list.txt`) containing the full paths to your VCF files, with one path per line. For example:
    
    ```text
    /path/to/sample1.vcf.gz
    /path/to/sample2.vcf
    /different/path/sample3.vcf.gz
    ```
    
    The sample names will be extracted from the filenames (without directories or extensions).

## Running the Workflow

### On a SLURM Cluster (Recommended for multiple samples)

1.  Navigate to the `snakemake` directory.
2.  Customize `run_variantcentrifuge.sh` if needed (e.g., SLURM account, partition).
3.  Submit the job:
    ```bash
    sbatch run_variantcentrifuge.sh variantcentrifuge.smk my_vc_config.yaml 20 path/to/slurm/profile
    ```
    Replace:
    -   `my_vc_config.yaml` with your adapted configuration file.
    -   `20` with the maximum number of parallel jobs you want Snakemake to manage.
    -   `path/to/slurm/profile` with the path to your Snakemake SLURM profile (optional).

SLURM logs for the main Snakemake process will be written to `slurm_vc_logs/`. Individual rule logs (from `variantcentrifuge` itself) will be in `base_output_folder/logs_vc/`.

### Locally (for testing or small datasets)

1.  Navigate to the `snakemake` directory.
2.  Activate your conda environment if `variantcentrifuge` is installed there.
3.  Run Snakemake directly:
    ```bash
    snakemake -s variantcentrifuge.smk --configfile my_vc_config.yaml --cores 8 --use-conda --conda-prefix path/to/conda_envs_dir
    ```
    Replace:
    -   `my_vc_config.yaml` with your config file.
    -   `8` with the number of cores to use locally.
    -   `path/to/conda_envs_dir` with a directory where Snakemake can store conda environments.

## Output Structure

The workflow will create the following structure within your specified `base_output_folder`:

```
base_output_folder/
├── SAMPLE_A/
│   ├── SAMPLE_A.vc_analysis.tsv  (Main TSV output)
│   ├── SAMPLE_A.vc_analysis.xlsx (Optional XLSX output)
│   ├── report/
│   │   └── index.html            (Optional HTML report)
│   │   └── variants.json
│   │   └── summary.json
│   │   └── igv/                  (Optional IGV reports)
│   │       └── ...
│   └── intermediate/             (Intermediate files from variantcentrifuge)
│       └── ...
│   └── bed_cache/                (BED cache from variantcentrifuge)
│       └── ...
├── SAMPLE_B/
│   └── ... (similar structure)
└── logs_vc/
    ├── SAMPLE_A.variantcentrifuge.log
    ├── SAMPLE_B.variantcentrifuge.log
    └── ...
```
The exact contents of the `SAMPLE_X` subdirectories depend on the flags passed to `variantcentrifuge` (e.g., `--xlsx`, `--html-report`).

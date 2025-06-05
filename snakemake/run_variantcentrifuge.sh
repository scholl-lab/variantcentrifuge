#!/bin/bash
#
#SBATCH --job-name=smk_vc_main
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=168:00:00 # Max walltime for the entire workflow
#SBATCH --mem=4000M      # Memory for the main Snakemake process
#SBATCH --output=slurm_vc_logs/%x-%j.log

################################################################################
# Usage:
#   sbatch run_variantcentrifuge.sh [SNAKEMAKE_FILE] [CONFIG_VC_FILE] [MAX_JOBS] [SLURM_PROFILE]
#
#   - SNAKEMAKE_FILE:   Path to the Snakemake workflow file
#                       (default: "variantcentrifuge.smk")
#   - CONFIG_VC_FILE:   Path to the VariantCentrifuge Snakemake config file 
#                       (default: "config_vc.yaml")
#   - MAX_JOBS:         Number of Snakemake jobs (in parallel) to use (default: 20)
#   - SLURM_PROFILE:    Path to a Snakemake SLURM profile (optional, e.g., "cubi-v1")
#                       If not provided, Snakemake runs locally or uses default cluster config.
################################################################################

# --- 1) Parse command-line arguments with defaults ---
SNAKEMAKE_FILE=${1:-"variantcentrifuge.smk"}
CONFIG_VC_FILE=${2:-"config_vc.yaml"}
MAX_JOBS=${3:-20}
SLURM_PROFILE_ARG=""
if [ -n "$4" ]; then
    SLURM_PROFILE_ARG="--profile=$4"
fi

# --- 2) HPC environment setup ---
# Setup TMPDIR for scratch space if not already set by scheduler
export TMPDIR=${TMPDIR:-"$HOME/scratch/tmp_vc/$$"} 
mkdir -p "$TMPDIR"
trap "rm -rf $TMPDIR" EXIT # Cleanup trap

# Create directory for SLURM logs from the main Snakemake process
mkdir -p slurm_vc_logs

echo "----------------------------------------------------"
echo "Starting VariantCentrifuge Snakemake Workflow"
echo "Timestamp: $(date)"
echo "Snakefile: $SNAKEMAKE_FILE"
echo "Config file: $CONFIG_VC_FILE"
echo "Max parallel jobs: $MAX_JOBS"
echo "SLURM Profile Arg: $SLURM_PROFILE_ARG"
echo "Temporary Directory: $TMPDIR"
echo "----------------------------------------------------"

# --- 3) Run Snakemake ---
# Ensure the script is run from the directory containing the Snakefile and config,
# or adjust paths accordingly. This assumes it's run from `snakemake/` directory.

snakemake \
    -s "$SNAKEMAKE_FILE" \
    --configfile "$CONFIG_VC_FILE" \
    -j "$MAX_JOBS" \
    --use-conda \
    --conda-prefix "$HOME/snakemake_conda_envs/vc_conda" \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    $SLURM_PROFILE_ARG \
    all # Target rule

EXIT_CODE=$?
echo "----------------------------------------------------"
echo "Snakemake workflow finished at $(date) with exit code $EXIT_CODE."
echo "----------------------------------------------------"
exit $EXIT_CODE

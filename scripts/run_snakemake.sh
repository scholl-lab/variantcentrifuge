#!/bin/bash
#SBATCH --job-name=sm_variantcentrifuge
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm_logs/%x-%j.log

# Auto-detecting launcher for the variantcentrifuge Snakemake workflow.
#
# Usage:
#   sbatch scripts/run_snakemake.sh [CONFIG_FILE] [EXTRA_SNAKEMAKE_ARGS...]
#   bash  scripts/run_snakemake.sh  [CONFIG_FILE] [EXTRA_SNAKEMAKE_ARGS...]
#
# The cluster environment (BIH / Charite / local) is detected automatically
# and the matching Snakemake profile from profiles/ is applied.

set -euo pipefail

# ── Cluster auto-detection ───────────────────────────────────────────────────
detect_cluster() {
    local fqdn
    fqdn=$(hostname -f 2>/dev/null || hostname)
    if [[ -d "/etc/xdg/snakemake/cubi-v1" ]] || [[ "${fqdn}" =~ cubi|bihealth ]]; then
        echo "bih"
    elif [[ -f "/etc/profile.d/conda.sh" ]] || [[ "${fqdn}" =~ charite|\.sc- ]]; then
        echo "charite"
    else
        echo "local"
    fi
}

CLUSTER=$(detect_cluster)
CONFIG_FILE="${1:-config/config_vc.yaml}"
shift 1 2>/dev/null || true

# ── Conda activation (Charite requires explicit sourcing) ────────────────────
if [[ "${CLUSTER}" == "charite" ]] && [[ -f /etc/profile.d/conda.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/conda.sh
fi
conda activate snakemake 2>/dev/null || true

# ── TMPDIR setup ─────────────────────────────────────────────────────────────
# Create a Snakemake-specific temp directory without overwriting TMPDIR if
# already set (e.g. by SLURM).
if [[ "${CLUSTER}" == "bih" ]]; then
    BASE_TMPDIR="${HOME}/scratch/tmp"
else
    BASE_TMPDIR="${TMPDIR:-/tmp}/snakemake"
fi
mkdir -p "${BASE_TMPDIR}"
SNAKEMAKE_TMPDIR=$(mktemp -d "${BASE_TMPDIR}/sm.XXXXXX")
export SNAKEMAKE_TMPDIR
if [[ -z "${TMPDIR:-}" ]]; then
    TMPDIR="${SNAKEMAKE_TMPDIR}"
    export TMPDIR
fi
trap 'rm -rf "${SNAKEMAKE_TMPDIR}"' EXIT

mkdir -p slurm_logs

# ── Launch ───────────────────────────────────────────────────────────────────
echo "=== variantcentrifuge Snakemake Launch ==="
echo "  Cluster:           ${CLUSTER}"
echo "  Config:            ${CONFIG_FILE}"
echo "  Profile:           profiles/${CLUSTER}"
echo "  TMPDIR:            ${TMPDIR:-<unset>}"
echo "  SNAKEMAKE_TMPDIR:  ${SNAKEMAKE_TMPDIR}"
echo "  Extra args:        $*"
echo "  Start:             $(date)"
echo "==========================================="

snakemake \
    -s workflow/Snakefile \
    --configfile "${CONFIG_FILE}" \
    --workflow-profile profiles/default \
    --profile "profiles/${CLUSTER}" \
    "$@"

echo "=== Finished: $(date) ==="

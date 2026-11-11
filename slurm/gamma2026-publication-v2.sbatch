#!/bin/bash
#SBATCH --job-name=ffr-gamma-gallery-v2
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/gamma2026-publication-v2-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
export MPLCONFIGDIR="$FFR_ROOT/rendered/mpl-gamma-v2-${SLURM_JOB_ID}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
mkdir -p "$MPLCONFIGDIR"
tools/gallery-env/bin/python -u scripts/gamma2026_publication_v2/publish.py "$@"

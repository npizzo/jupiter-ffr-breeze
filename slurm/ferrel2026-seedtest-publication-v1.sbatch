#!/bin/bash
#SBATCH --job-name=ffr-seed-publish-v1
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ferrel2026-seedtest-publication-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
: "${SLURM_JOB_ID:?Reader and rendering require an allocated Unity job}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
mkdir -p "$MPLCONFIGDIR"
tools/gallery-env/bin/python -u scripts/ferrel2026_seedtest_publication_v1/report.py "$@"

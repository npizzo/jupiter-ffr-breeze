#!/bin/bash
#SBATCH --job-name=ffr-ferrel-publish
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ferrel2026-publication-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
mkdir -p "$MPLCONFIGDIR"
tools/gallery-env/bin/python -u scripts/ferrel2026_publication/render_report.py "$@"

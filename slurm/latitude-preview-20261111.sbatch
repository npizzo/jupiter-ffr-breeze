#!/bin/bash
#SBATCH --job-name=ffr-latitude-preview
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/latitude-preview-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/latitude_preview_20261111/prepare.py "$1"

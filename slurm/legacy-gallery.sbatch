#!/bin/bash
#SBATCH --job-name=ffr-legacy-gallery
#SBATCH --partition=uri-cpu,cpu
#SBATCH --array=0-3
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/legacy-gallery-%A_%a.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/render_legacy_archive.py

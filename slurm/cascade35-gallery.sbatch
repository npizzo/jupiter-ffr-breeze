#!/bin/bash
#SBATCH --job-name=ffr-beta35-gallery
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cascade35-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/cascade35/prepare_gallery.py "$1"

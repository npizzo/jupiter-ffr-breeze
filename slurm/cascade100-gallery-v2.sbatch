#!/bin/bash
#SBATCH --job-name=ffr-gallery100
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --output=logs/cascade100-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/prepare_cascade100_gallery_v2.py

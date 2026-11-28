#!/bin/bash
#SBATCH --job-name=nov28-gallery
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/nov28-gallery-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-nov28-$1"
tools/point4-analysis-env/bin/python scripts/nov28_completed_publication_v1/render.py "$1"

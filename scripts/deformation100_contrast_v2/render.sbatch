#!/bin/bash
#SBATCH --job-name=deformation100-contrast
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/deformation100-contrast-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-deformation100-contrast"
tools/point4-analysis-env/bin/python scripts/deformation100_contrast_v2/render.py

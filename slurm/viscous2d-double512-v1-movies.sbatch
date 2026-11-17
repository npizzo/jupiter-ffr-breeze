#!/bin/bash
#SBATCH --job-name=ffr-512-double-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/viscous2d-double512-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/viscous2d_double512_v1/render.py data/viscous2d-long512-v1-1848831 data/viscous2d-fast512-v1-1872263 "$1" "$2"

#!/bin/bash
#SBATCH --job-name=ffr-double160-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/viscous2d-double512-restart-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/viscous2d_double512_restart_v1/render160.py data/viscous2d-long512-v1-1848831 data/viscous2d-fast512-v1-1872263 "$1" "$2"

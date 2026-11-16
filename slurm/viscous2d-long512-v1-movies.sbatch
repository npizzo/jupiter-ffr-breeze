#!/bin/bash
#SBATCH --job-name=ffr-512-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/viscous2d-long512-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/viscous2d_long512_v1/render.py "$1" "$2"

#!/bin/bash
#SBATCH --job-name=ffr-feedback-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/viscous2d-feedback512-movies-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/viscous2d_feedback512_v1/render.py "$1" "$2"

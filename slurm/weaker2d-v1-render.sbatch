#!/bin/bash
#SBATCH --job-name=ffr-2d-weak-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/weaker2d-v1-render-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/weaker2d_v1/analyze_render.py data/instability2d-long-v1-1817684 "$1" "$2"

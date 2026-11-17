#!/bin/bash
#SBATCH --job-name=ffr-stress-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/stress2d-movies-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/stress2d_v1/render.py "$1" "$2"

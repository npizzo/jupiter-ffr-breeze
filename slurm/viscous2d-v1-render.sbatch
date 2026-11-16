#!/bin/bash
#SBATCH --job-name=ffr-2d-visc-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/viscous2d-v1-render-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/viscous2d_v1/render.py "$1" "$2"

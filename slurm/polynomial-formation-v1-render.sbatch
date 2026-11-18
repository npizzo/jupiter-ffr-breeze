#!/bin/bash
#SBATCH --job-name=ffr-poly-form-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/polynomial-formation-render-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/polynomial_formation_v1/render.py "$1" "$2"

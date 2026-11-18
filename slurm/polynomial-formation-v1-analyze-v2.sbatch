#!/bin/bash
#SBATCH --job-name=ffr-poly-form-analysis
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/polynomial-formation-analyze-%j.out
#SBATCH --no-requeue
set -euo pipefail
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/polynomial_formation_v1/analyze_v2.py data/viscous2d-feedback512-v1-1925368 "$1" "$2"

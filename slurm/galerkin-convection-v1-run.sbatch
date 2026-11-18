#!/bin/bash
#SBATCH --job-name=ffr-galerkin-nonlinear
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/galerkin-convection-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$1/validation/PASS.json"
tools/point4-analysis-env/bin/python scripts/galerkin_convection_v1/campaign.py "$1"

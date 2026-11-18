#!/bin/bash
#SBATCH --job-name=ffr-galerkin-render
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=logs/galerkin-convection-render-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
test -f "$1/CAMPAIGN_COMPLETE.json"
tools/point4-analysis-env/bin/python scripts/galerkin_convection_v1/analyze.py "$1" "$2"

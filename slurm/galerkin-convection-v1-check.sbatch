#!/bin/bash
#SBATCH --job-name=ffr-galerkin-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/galerkin-convection-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/galerkin_convection_v1/check_campaign.py "data/galerkin-convection-v1-${SLURM_JOB_ID}"

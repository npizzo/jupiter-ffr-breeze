#!/bin/bash
#SBATCH --job-name=ffr-mode-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/convective-response-modes-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/convective_response_v1/optimized_screen.py experiments/convective-response-v1/reference-profile.json "$1/report.json" "data/convective-response-modes-${SLURM_JOB_ID}"

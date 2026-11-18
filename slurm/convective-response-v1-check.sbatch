#!/bin/bash
#SBATCH --job-name=ffr-response-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/convective-response-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
test -f "$1/PASS.json"
tools/point4-analysis-env/bin/python scripts/convective_response_v1/check_response.py "$1/basis.npz" "data/convective-response-check-${SLURM_JOB_ID}"

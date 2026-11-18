#!/bin/bash
#SBATCH --job-name=ffr-response
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=logs/convective-response-run-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/convective_response_v1/run.py "$1/basis.npz" "$3" --validation "$2/PASS.json"

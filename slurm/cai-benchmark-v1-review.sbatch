#!/bin/bash
#SBATCH --job-name=ffr-cai-review-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=logs/cai-benchmark-v1-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec tools/point4-analysis-env/bin/python scripts/cai_benchmark_v1/review.py "rendered/cai-benchmark-v1-$SLURM_JOB_ID" "$@"

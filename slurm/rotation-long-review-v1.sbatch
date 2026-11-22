#!/bin/bash
#SBATCH --job-name=ffr-rotation-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/rotation-long-review-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="/tmp/ffr-mpl-$SLURM_JOB_ID"
exec tools/point4-analysis-env/bin/python scripts/rotation_long_review_v1/review.py "data/rotation-long-review-$SLURM_JOB_ID/publication"

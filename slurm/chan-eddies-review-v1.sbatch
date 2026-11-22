#!/bin/bash
#SBATCH --job-name=ffr-eddies-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/chan-eddies-review-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="/tmp/ffr-mpl-$SLURM_JOB_ID"
exec tools/point4-analysis-env/bin/python scripts/chan_eddies_review_v1/review.py "data/chan-eddies-review-$SLURM_JOB_ID/publication"

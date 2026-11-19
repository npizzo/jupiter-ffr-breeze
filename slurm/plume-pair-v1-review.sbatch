#!/bin/bash
#SBATCH --job-name=ffr-plume-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/plume-pair-v1-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/plume-mpl-$SLURM_JOB_ID
tools/point4-analysis-env/bin/python scripts/plume_pair_v1/review.py "rendered/plume-pair-v1-$SLURM_JOB_ID" "$1" "$2"

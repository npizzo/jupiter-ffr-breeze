#!/bin/bash
#SBATCH --job-name=ffr-sign-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fplane-sign-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/fplane_sign_v1/review.py data/fplane-sign-v1 rendered/fplane-sign-v1

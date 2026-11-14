#!/bin/bash
#SBATCH --job-name=ffr-resolution-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ferrel35-resolution-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python scripts/ferrel35_resolution_review_v1/review.py "$1" "$2"
tools/gallery-env/bin/python scripts/ferrel35_resolution_review_v1/render.py "$1" "$2"

#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-face-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/ferrel35-face-review-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python scripts/ferrel35_face_v2/review.py "$1/production" "$2"
tools/gallery-env/bin/python scripts/ferrel35_face_v2/render.py "$1/production" "$2"

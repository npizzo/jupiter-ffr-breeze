#!/bin/bash
#SBATCH --job-name=ffr-wide100-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ferrel35-wide100-movie-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python scripts/ferrel35_width100_gallery/prepare.py wide data/ferrel35-width-v4 rendered/ferrel35-wide100-v1

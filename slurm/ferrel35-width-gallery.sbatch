#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-width-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ferrel35-width-movie-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python scripts/ferrel35_width_gallery/prepare.py "$1" "$2" "$3"

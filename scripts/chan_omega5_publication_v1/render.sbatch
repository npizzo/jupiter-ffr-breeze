#!/bin/bash
#SBATCH --job-name=chan-omega5-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/chan-omega5-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-nov26-$1"
tools/point4-analysis-env/bin/python scripts/chan_omega5_publication_v1/render.py "$1"

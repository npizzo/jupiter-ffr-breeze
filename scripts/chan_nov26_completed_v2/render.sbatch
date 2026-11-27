#!/bin/bash
#SBATCH --job-name=chan-completed-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=logs/chan-nov26-completed-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-nov26-$1"
tools/point4-analysis-env/bin/python scripts/chan_nov26_completed_v2/render.py "$1"

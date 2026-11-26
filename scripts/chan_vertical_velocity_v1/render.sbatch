#!/bin/bash
#SBATCH --job-name=chan-vertical-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/chan-vertical-movies-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-chan-vertical-v1"
tools/point4-analysis-env/bin/python scripts/chan_vertical_velocity_v1/render.py

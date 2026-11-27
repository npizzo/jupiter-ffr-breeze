#!/bin/bash
#SBATCH --job-name=chan-tracer-contrast-v3
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-tracer-contrast-v3-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-tracer-v3-$1"
tools/point4-analysis-env/bin/python scripts/chan_tracer_contrast_v3/render.py "$1"

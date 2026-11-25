#!/bin/bash
#SBATCH --job-name=chan-totalE2000-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-total-energy2000-render-v1-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-totalE2000"
tools/point4-analysis-env/bin/python scripts/chan_total_energy2000_publication_v1/render.py

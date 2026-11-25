#!/bin/bash
#SBATCH --job-name=chan-anel6000-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=logs/chan-clean-anelastic6000-render-v1-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-anel6000"
tools/point4-analysis-env/bin/python scripts/chan_clean_anelastic6000_publication_v1/render.py

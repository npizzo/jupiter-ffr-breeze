#!/bin/bash
#SBATCH --job-name=chan-anel4000-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/chan-clean-anelastic4000-render-v1-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-anel4000"
tools/point4-analysis-env/bin/python scripts/chan_clean_anelastic4000_publication_v1/render.py

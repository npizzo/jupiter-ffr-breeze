#!/bin/bash
#SBATCH --job-name=cai-anelastic-segment-movie
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/cai-anelastic-segment-movie-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR="$PWD/work/mpl-cai-anelastic-segment"
tools/point4-analysis-env/bin/python scripts/cai_anelastic_segment_publication_v1/render.py

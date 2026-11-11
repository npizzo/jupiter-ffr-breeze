#!/bin/bash
#SBATCH --job-name=ffr-earlier-gallery
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/earlier-gallery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/render_earlier_archive.py

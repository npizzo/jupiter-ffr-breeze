#!/bin/bash
#SBATCH --job-name=ffr-archive-web
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/earlier-gallery-web-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/finalize_earlier_archive.py

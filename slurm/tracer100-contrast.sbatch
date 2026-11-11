#!/bin/bash
#SBATCH --job-name=ffr-tracer100-contrast
#SBATCH --partition=uri-cpu,cpu
#SBATCH --array=0-2
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/tracer100-contrast-%A_%a.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
export MPLCONFIGDIR="$PWD/rendered/mplconfig"
tools/gallery-env/bin/python -u scripts/tracer100_contrast/prepare.py

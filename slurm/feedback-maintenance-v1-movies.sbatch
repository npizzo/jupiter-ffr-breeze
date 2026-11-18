#!/bin/bash
#SBATCH --job-name=ffr-maint-movies
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/feedback-maint-movies-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/feedback_maintenance_v1/render.py "$1" "$2"

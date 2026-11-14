#!/bin/bash
#SBATCH --job-name=ffr-point4-analysis-gate
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fplane-point4-recovery-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
 tools/point4-analysis-env/bin/python scripts/fplane_point4_recovery_v3/validate.py

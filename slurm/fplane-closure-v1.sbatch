#!/bin/bash
#SBATCH --job-name=ffr-convective-closure
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fplane-closure-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
tools/point4-analysis-env/bin/python scripts/fplane_closure_v1/analyze.py data "rendered/fplane-closure-v1-${SLURM_JOB_ID}"

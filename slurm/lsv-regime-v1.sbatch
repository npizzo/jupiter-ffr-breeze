#!/bin/bash
#SBATCH --job-name=ffr-lsv-regime-audit
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/lsv-regime-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec tools/point4-analysis-env/bin/python scripts/lsv_regime_v1/audit.py "rendered/lsv-regime-v1-$SLURM_JOB_ID"

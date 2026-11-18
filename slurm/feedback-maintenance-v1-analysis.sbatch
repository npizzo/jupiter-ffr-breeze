#!/bin/bash
#SBATCH --job-name=ffr-maint-analysis
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/feedback-maint-analysis-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python_bin=tools/point4-analysis-env/bin/python
"$python_bin" scripts/feedback_maintenance_v1/analyze_existing.py data/viscous2d-feedback512-v1-1925368 "rendered/feedback-maintenance-existing-${SLURM_JOB_ID}"

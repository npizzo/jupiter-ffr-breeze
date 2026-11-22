#!/bin/bash
#SBATCH --job-name=ffr-anelastic-results
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/omega3-anelastic-results-v1-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4
tools/point4-analysis-env/bin/python scripts/omega3_anelastic_results_v1/review.py "data/omega3-anelastic-results-v1-$SLURM_JOB_ID"

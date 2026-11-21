#!/bin/bash
#SBATCH --job-name=ffr-anelastic6000-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:15:00
#SBATCH --output=logs/anelastic-continuation-review-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MPLCONFIGDIR="/tmp/ffr-mpl-$SLURM_JOB_ID"
exec tools/point4-analysis-env/bin/python scripts/anelastic_continuation_review_v1/review.py data/anelastic-sequence-v3-launchfix-2058241 data/anelastic-sequence-v3-extendfix-2071578 "data/anelastic-continuation-review-$SLURM_JOB_ID/publication"

#!/bin/bash
#SBATCH --job-name=ffr-anelastic12000-review
#SBATCH --partition=uri-cpu,cpu
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/anelastic-review12000-v1-%j.out
#SBATCH --no-requeue
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=4 MPLCONFIGDIR="/tmp/ffr-mpl-$SLURM_JOB_ID"
exec tools/point4-analysis-env/bin/python scripts/anelastic_review12000_v1/review.py data/anelastic-sequence-v3-launchfix-2058241 data/anelastic-sequence-v3-extendfix-2071578 data/anelastic-sequence-v3-to12000-2077087 "data/anelastic-review12000-$SLURM_JOB_ID/publication"

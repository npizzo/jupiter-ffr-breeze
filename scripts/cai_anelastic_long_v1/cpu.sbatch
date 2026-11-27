#!/bin/bash
#SBATCH --job-name=cai-anelastic-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:20:00
#SBATCH --output=logs/cai-anelastic-cpu-%j.out
set -euo pipefail
source local.env
source scripts/cai_anelastic_long_v1/environment.sh
python scripts/cai_anelastic_long_v1/gate.py preflight "runs/cai-anelastic-long-v1/cpu-$SLURM_JOB_ID"

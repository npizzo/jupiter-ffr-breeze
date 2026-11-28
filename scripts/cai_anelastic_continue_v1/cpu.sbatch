#!/bin/bash
#SBATCH --job-name=cai-anelastic-resume-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=00:20:00
#SBATCH --output=logs/cai-anelastic-resume-cpu-%j.out
set -euo pipefail
source local.env
source scripts/cai_anelastic_continue_v1/environment.sh
python scripts/cai_anelastic_continue_v1/gate.py preflight "runs/cai-anelastic-continue-v1/cpu-$SLURM_JOB_ID"

#!/bin/bash
#SBATCH --job-name=marcus-fourday-cpu
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:25:00
#SBATCH --output=logs/marcus-fourday-cpu-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_fourday_v1/environment.sh
python scripts/marcus_cv1_fourday_v1/gate.py preflight "runs/marcus-cv1-fourday-v1/cpu-$SLURM_JOB_ID"

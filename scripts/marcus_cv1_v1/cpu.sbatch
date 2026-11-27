#!/bin/bash
#SBATCH --job-name=marcus-equation-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:25:00
#SBATCH --output=logs/marcus-cv1-cpu-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_v1/environment.sh
python scripts/marcus_cv1_v1/gate.py preflight "runs/marcus-cv1-v1/cpu-$SLURM_JOB_ID"

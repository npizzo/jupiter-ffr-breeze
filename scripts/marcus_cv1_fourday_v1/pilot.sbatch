#!/bin/bash
#SBATCH --job-name=marcus-resolved-pilot
#SBATCH --partition=uri-gpu,gpu-highvram,gpu-preempt-highvram
#SBATCH --constraint=(a100|h100)&vram80
#SBATCH --exclude=uri-gpu006
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --no-requeue
#SBATCH --output=logs/marcus-fourday-pilot-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_fourday_v1/environment.sh
python scripts/marcus_cv1_fourday_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/marcus_cv1_fourday_v1/driver.jl pilot "runs/marcus-cv1-fourday-v1/pilot-$SLURM_JOB_ID"

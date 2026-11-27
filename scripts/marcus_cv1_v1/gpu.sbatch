#!/bin/bash
#SBATCH --job-name=marcus-fixture-gpu
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/marcus-cv1-gpu-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_v1/environment.sh
python scripts/marcus_cv1_v1/gate.py verify "$1"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/marcus_cv1_v1/driver.jl gpu_check "runs/marcus-cv1-v1/gpu-$SLURM_JOB_ID"

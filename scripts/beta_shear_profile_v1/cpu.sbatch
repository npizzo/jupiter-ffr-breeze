#!/bin/bash
#SBATCH --job-name=ffr-repair-preflight
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/beta-shear-profile-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
out="data/beta-shear-profile-v1-cpu-$SLURM_JOB_ID"
python scripts/beta_shear_profile_v1/gate.py preflight "$out"
sbatch --parsable scripts/beta_shear_profile_v1/gpu.sbatch "$out/LAUNCH_PASS.json" run

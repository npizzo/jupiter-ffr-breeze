#!/bin/bash
#SBATCH --job-name=ffr-shear-physics-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:25:00
#SBATCH --output=logs/beta-shear-physics-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
out="data/beta-shear-physics-v1-cpu-$SLURM_JOB_ID"
python scripts/beta_shear_physics_v1/gate.py preflight "$out"
python scripts/beta_shear_physics_v1/gate.py verify "$out/LAUNCH_PASS.json"
sbatch --parsable scripts/beta_shear_physics_v1/gpu.sbatch "$out/LAUNCH_PASS.json" complete

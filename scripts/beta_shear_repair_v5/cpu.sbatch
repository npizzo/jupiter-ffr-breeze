#!/bin/bash
#SBATCH --job-name=ffr-repair-preflight
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/beta-shear-repair-v5-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
out="data/beta-shear-repair-v5-cpu-$SLURM_JOB_ID"
python scripts/beta_shear_repair_v5/gate.py preflight "$out"
sbatch --parsable scripts/beta_shear_repair_v5/gpu.sbatch "$out/LAUNCH_PASS.json" run
sbatch --parsable scripts/beta_shear_repair_v5/gpu.sbatch "$out/LAUNCH_PASS.json" forensic

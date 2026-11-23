#!/bin/bash
#SBATCH --job-name=ffr-jet-diag-preflight
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:10:00
#SBATCH --output=logs/beta-shear-diag-v2-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
out="data/beta-shear-diag-v2-cpu-$SLURM_JOB_ID"
python scripts/beta_shear_diagnosis_v2/gate.py preflight "$out"
sbatch --parsable scripts/beta_shear_diagnosis_v2/gpu.sbatch "$out/LAUNCH_PASS.json"

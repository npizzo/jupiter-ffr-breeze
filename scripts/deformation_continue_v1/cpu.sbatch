#!/bin/bash
#SBATCH --job-name=ffr-deformation-restart-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:20:00
#SBATCH --output=logs/deformation-continue-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/deformation_continue_v1/gate.py preflight "data/deformation-continue-v1-cpu-$SLURM_JOB_ID"
python scripts/deformation_continue_v1/release.py "data/deformation-continue-v1-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"

#!/bin/bash
#SBATCH --job-name=ffr-deformation-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH --output=logs/deformation-pair-v1-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/deformation_pair_v1/gate.py preflight "data/deformation-pair-v1-cpu-$SLURM_JOB_ID"
python scripts/deformation_pair_v1/release.py "data/deformation-pair-v1-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"

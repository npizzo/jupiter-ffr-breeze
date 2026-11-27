#!/bin/bash
#SBATCH --job-name=ffr-convective-cap-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/convective-deformation-v2-cpu-%j.out
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
python scripts/convective_deformation_v2/gate.py preflight "data/convective-deformation-v2-cpu-$SLURM_JOB_ID"
python scripts/convective_deformation_v2/release.py "data/convective-deformation-v2-cpu-$SLURM_JOB_ID/LAUNCH_PASS.json"

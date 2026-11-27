#!/bin/bash
#SBATCH --job-name=cai-anelastic-B7
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=a100&vram80
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cai-anelastic-long-%j.out
set -euo pipefail
source local.env
source scripts/cai_anelastic_long_v1/environment.sh
python scripts/cai_anelastic_long_v1/gate.py verify "$1"
data_root=/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127
out="$data_root/anelastic-B7-$SLURM_JOB_ID"
mkdir -p runs/cai-anelastic-long-v1
printf '%s\n' "$out" > "runs/cai-anelastic-long-v1/job-$SLURM_JOB_ID-output.txt"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/cai_anelastic_long_v1/driver.jl run "$out"

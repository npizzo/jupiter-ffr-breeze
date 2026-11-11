#!/bin/bash
#SBATCH --job-name=ffr-seedv1-1024
#SBATCH --partition=uri-gpu
#SBATCH --gpus=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel2026-seedtest-v1-1024-%j.out
#SBATCH --no-requeue
set -euo pipefail
test "$#" -eq 4
test "$1" = modes1024_b162
bash scripts/ferrel2026_seedtest_v1/gpu_job.sh "$4" scripts/ferrel2026_seedtest_v1/campaign.jl "$@"

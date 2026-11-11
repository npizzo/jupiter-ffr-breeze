#!/bin/bash
#SBATCH --job-name=ffr-seedv1-512
#SBATCH --partition=uri-gpu
#SBATCH --gpus=l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ferrel2026-seedtest-v1-512-%j.out
#SBATCH --no-requeue
set -euo pipefail
test "$#" -eq 4
case_name="$1"
case "$case_name" in zero512_b0|zero512_b162|modes512_b162) ;; *) exit 2 ;; esac
bash scripts/ferrel2026_seedtest_v1/gpu_job.sh "$4" scripts/ferrel2026_seedtest_v1/campaign.jl "$@"

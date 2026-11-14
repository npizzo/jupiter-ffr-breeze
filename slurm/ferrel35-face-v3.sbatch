#!/bin/bash
#SBATCH --job-name=ffr-ferrel35-face300
#SBATCH --partition=uri-gpu
#SBATCH --gpus=1
#SBATCH --constraint=vram48
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=logs/ferrel35-face-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
"$JULIA_BIN" --project=. scripts/ferrel35_face_v3/cuda_launch.jl "$1"

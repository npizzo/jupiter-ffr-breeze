#!/bin/bash
#SBATCH --job-name=cai-B-benchmark
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cai-benchmark-%j.out
set -euo pipefail
source ./local.env
source scripts/cai2022_benchmark_v2/environment.sh
python scripts/cai2022_benchmark_v2/gate.py verify "$1"
case_name="$2"
data_root=/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/cai2022_benchmark_v2/driver.jl gpu_validate "$case_name" "$data_root/gpu_validate-$case_name-$SLURM_JOB_ID"

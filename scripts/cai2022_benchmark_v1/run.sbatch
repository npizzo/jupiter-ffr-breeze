#!/bin/bash
#SBATCH --job-name=cai-B-benchmark
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=fp64
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cai-benchmark-%j.out
set -euo pipefail
source ./local.env
source scripts/cai2022_benchmark_v1/environment.sh
python scripts/cai2022_benchmark_v1/gate.py verify "$1"
case_name="$2"
data_root=/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/cai2022_benchmark_v1/driver.jl run "$case_name" "$data_root/run-$case_name-$SLURM_JOB_ID"

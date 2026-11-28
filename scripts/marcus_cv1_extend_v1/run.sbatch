#!/bin/bash
#SBATCH --job-name=marcus-CV1-extend
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=(a100|h100)&vram80
#SBATCH --exclude=uri-gpu006
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --time-min=04:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/marcus-extend-%j.out
set -euo pipefail
source local.env
source scripts/marcus_cv1_extend_v1/environment.sh
out="/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127/marcus-CV1-extend-$SLURM_JOB_ID"
case_name="$3"
on_error() {
    code=$?
    python scripts/marcus_cv1_extend_v1/chain.py failed "$case_name" "$out" "$SLURM_JOB_ID" "$code"
    exit "$code"
}
trap on_error ERR
python scripts/marcus_cv1_extend_v1/gate.py verify "$1"
python scripts/marcus_cv1_extend_v1/chain.py verify-request "$3" "$2"
python scripts/marcus_cv1_extend_v1/chain.py running "$3" "$out" "$SLURM_JOB_ID"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/marcus_cv1_extend_v1/driver.jl resume "$out" "$2" "$3"
python scripts/marcus_cv1_extend_v1/chain.py finish "$3" "$out" "$1"

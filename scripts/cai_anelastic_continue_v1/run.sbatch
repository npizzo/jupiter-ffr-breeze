#!/bin/bash
#SBATCH --job-name=cai-B7-continue5000
#SBATCH --partition=uri-gpu,gpu-highvram
#SBATCH --constraint=(a100|h100)&vram80
#SBATCH --exclude=uri-gpu006
#SBATCH --gpus=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --time-min=01:00:00
#SBATCH --no-requeue
#SBATCH --output=logs/cai-anelastic-continue-%j.out
set -euo pipefail
source local.env
source scripts/cai_anelastic_continue_v1/environment.sh
out="/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127/anelastic-B7-cont-$SLURM_JOB_ID"
on_error() {
    code=$?
    python scripts/cai_anelastic_continue_v1/chain.py failed "$out" "$SLURM_JOB_ID" "$code"
    exit "$code"
}
trap on_error ERR
python scripts/cai_anelastic_continue_v1/gate.py verify "$1"
python scripts/cai_anelastic_continue_v1/chain.py verify-request "$2"
python scripts/cai_anelastic_continue_v1/chain.py running "$out" "$SLURM_JOB_ID"
"$JULIA_BIN" --pkgimages=no --project="$TE_PROJECT" scripts/cai_anelastic_continue_v1/driver.jl resume "$out" "$2"
python scripts/cai_anelastic_continue_v1/chain.py finish "$out" "$1"

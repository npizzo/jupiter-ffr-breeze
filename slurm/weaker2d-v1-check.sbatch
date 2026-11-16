#!/bin/bash
#SBATCH --job-name=ffr-2d-weak-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/weaker2d-v1-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source ./local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
out="data/weaker2d-v1-check-${SLURM_JOB_ID}"
mkdir "$out"
source_state=data/instability2d-light-v1-1810727/checkpoint.jld2
"$JULIA_BIN" --startup-file=no -e 'for p in ARGS; ex=Meta.parseall(read(p,String)); any(x -> x isa Expr && x.head in (:error,:incomplete), ex.args) && error("Parse failure in $p"); end; println("SOURCE_PARSE_PASS")' scripts/weaker2d_v1/{run,check}.jl
"$JULIA_BIN" --project=. scripts/weaker2d_v1/check.jl "$source_state" "$out/scaling.toml"
for case_name in forced random_sign cyclonic_sign; do
 "$JULIA_BIN" --project=. scripts/weaker2d_v1/run.jl "$out/$case_name" "$case_name" "$source_state" probe
done
tools/point4-analysis-env/bin/python scripts/weaker2d_v1/verify_probes.py "$out"

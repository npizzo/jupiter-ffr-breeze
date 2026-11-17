#!/bin/bash
#SBATCH --job-name=ffr-stress-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/stress2d-check-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --startup-file=no -e 'function check(ex); ex isa Expr || return; ex.head in (:error,:incomplete) && error(ex); foreach(check,ex.args); end; foreach(p->check(Meta.parseall(read(p,String))),ARGS); println("PARSE_PASS")' scripts/stress2d_v1/model.jl scripts/stress2d_v1/reference.jl scripts/stress2d_v1/run.jl scripts/stress2d_v1/check.jl
"$JULIA_BIN" --project=. scripts/stress2d_v1/check.jl data/viscous2d-refine-v1-1848240/initial.jld2 "data/stress2d-check-${SLURM_JOB_ID}"

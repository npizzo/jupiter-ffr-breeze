#!/bin/bash
#SBATCH --job-name=ffr-stress-driver-check
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/stress2d-release-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
CHECK=$1
"$JULIA_BIN" --project=. scripts/stress2d_v1/run.jl "$CHECK/driver-smoke" "$CHECK/initial.jld2" 1 "$CHECK/PASS.toml" biased_plus
tools/point4-analysis-env/bin/python - "$CHECK" <<'PY'
from pathlib import Path
import sys,tomllib,csv,json
p=Path(sys.argv[1])/'driver-smoke'
s=tomllib.loads((p/'status.toml').read_text());assert s['status']=='COMPLETED' and abs(s['t']-1)<1e-8
frames=list(csv.DictReader((p/'frames.csv').open()));assert [float(x['t']) for x in frames]==[0,1]
assert all((p/x['file']).stat().st_size==256*256*4 for x in frames)
assert (p/'checkpoint.jld2').is_file()
(p.parent/'DRIVER_PASS.json').write_text(json.dumps({'status':'PASS','frames':2,'end_time':1})+'\n')
print('DRIVER_PASS')
PY

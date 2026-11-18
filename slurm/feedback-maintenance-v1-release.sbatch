#!/bin/bash
#SBATCH --job-name=ffr-maint-driver
#SBATCH --partition=uri-cpu,cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/feedback-maint-driver-%j.out
#SBATCH --no-requeue
set -euo pipefail
source local.env
source scripts/cuda_repair/environment.sh
export JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
"$JULIA_BIN" --project=. scripts/feedback_maintenance_v1/run.jl "$1/driver-smoke" data/viscous2d-feedback512-v1-1925368/biased_plus/to_400/checkpoint.jld2 400.5 "$1/PASS.toml" no_low_quadratic
tools/point4-analysis-env/bin/python - "$1" <<'PY'
import sys,csv,tomllib,json
from pathlib import Path
p=Path(sys.argv[1]);o=p/'driver-smoke'
s=tomllib.loads((o/'status.toml').read_text());assert s['status']=='COMPLETED' and s['t']==400.5
fr=list(csv.DictReader((o/'frames.csv').open()));assert [float(r['t']) for r in fr]==[400.,400.5]
sc=list(csv.DictReader((o/'scalars.csv').open()));assert [float(r['t']) for r in sc]==[400.,400.25,400.5]
assert all(abs(float(r['quadratic_lowk_energy_rate']))<1e-12 for r in sc)
assert (o/'checkpoint.jld2').stat().st_size>0
(p/'DRIVER_PASS.json').write_text(json.dumps({'status':'PASS','target':400.5,'frames':2,'samples':3})+'\n')
print('DRIVER_PASS')
PY

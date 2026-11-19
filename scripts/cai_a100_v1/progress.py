"""Unity-only launcher/elapsed-time reporter. Reads logs, never model fields."""
import csv,json,os,signal,subprocess,sys,time,tomllib
from pathlib import Path
config,out=sys.argv[1:];root=Path(out)
cmd=[os.environ['JULIA_BIN'],'--project=.','scripts/cai_a100_v1/launch.jl',config,out]
started=float(os.environ['FFR_JOB_START_EPOCH'])
child=subprocess.Popen(cmd)
def forward(sig,frame):
 if child.poll() is None:child.send_signal(sig)
signal.signal(signal.SIGTERM,forward);signal.signal(signal.SIGINT,forward)
def record():
 path=root/'diagnostics.csv'
 if not path.exists():return
 try:
  with path.open() as h:rows=list(csv.DictReader(h))
  rows=[r for r in rows if r.get('wall_w_max') is not None]
  if not rows:return
  last=rows[-1];t=float(last['time']);elapsed=time.time()-started
  report=dict(case=Path(config).stem,simulated_time=t,target_time=5000,job_elapsed_seconds=elapsed,allocation_limit_seconds=14400,timing_is_informational=True,status='RUNNING' if child.poll() is None else 'PROCESS_EXITED',process_exit_code=child.poll())
  late=[r for r in rows if float(r['time'])>=100][-7:]
  if len(late)>=4:
   dt=float(late[-1]['time'])-float(late[0]['time']);dw=float(late[-1]['wall_seconds'])-float(late[0]['wall_seconds'])
   rate=dw/dt
   report.update(seconds_per_simulation_unit=rate,forecast_total_seconds=elapsed+max(0,5000-t)*rate,forecast_total_seconds_with_margin=elapsed+1.25*max(0,5000-t)*rate,forecast_note='Measured in this production allocation after t>=100; later CFL can change. Hard wall limit remains unchanged.')
  if (root/'metadata.toml').exists():
   meta=tomllib.loads((root/'metadata.toml').read_text());report['model_status']=meta['status']
  temp=root/'progress.tmp';temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(root/'progress.json')
 except (ValueError,KeyError,OSError,tomllib.TOMLDecodeError):
  # A concurrently written diagnostic line will be retried on the next report.
  return
while child.poll() is None:
 record();time.sleep(30)
record();sys.exit(child.returncode)

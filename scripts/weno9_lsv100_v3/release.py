"""Submit the four authorized t=10 -> 100 continuations only after exact-driver preflight."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys,tomllib
root=Path(__file__).resolve().parents[2]
here=Path(__file__).resolve().parent
folder=root/'runs/weno9-lsv100-v3';folder.mkdir(parents=True,exist_ok=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 record=folder/'submission.json'
 if record.exists():print(record.read_text());sys.exit(0)
 receipt=sys.argv[1]
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 for i,(eq,jet) in enumerate([('compressible',0.),('compressible',.015625),('anelastic',0.),('anelastic',.015625)]):
  p=root/f'data/weno9-startup-v2-pilot-2177387_{i}/parent'
  m=tomllib.loads((p/'metadata.toml').read_text())
  assert m['status']=='COMPLETE' and m['final_time']==10 and m['equations']==eq and m['jet_speed']==jet
  assert [m[k] for k in ['nx','ny','nz']]==[128,256,80]
  assert m['closure']=='none' and m['advection_order']==9 and m['initial_rms']==.136
  assert len(list((p/'checkpoints').glob('*.jld2')))==1
  for f in ['health_reference.jld2','budget_reference.jld2']:assert (p/f).is_file()
 assert shutil.disk_usage(root).free>20*1024**3
 job=subprocess.check_output(['sbatch','--parsable',str(here/'run.sbatch'),receipt,'2177387'],cwd=root,text=True).strip()
 state={'status':'SUBMITTED','gpu_array_job':job,'receipt':receipt,'restart_time':10,'target_time':100,'parent_array':2177387,'automatic_further_extension':False}
 record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

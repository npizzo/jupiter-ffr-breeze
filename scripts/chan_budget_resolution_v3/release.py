"""Release baseline first; unlock the other two only after its GPU checks pass."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys,tomllib
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
folder=root/'runs/chan-budget-resolution-v3';folder.mkdir(parents=True,exist_ok=True)
receipt=sys.argv[1];phase=sys.argv[2] if len(sys.argv)>2 else 'first'
assert phase in ('first','remaining')
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 record=folder/('submission-'+phase+'.json')
 if record.exists():print(record.read_text());sys.exit(0)
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 if phase=='remaining':
  check=Path(sys.argv[3])/'GPU_PASS.toml'
  assert tomllib.loads(check.read_text())['status']=='PASS'
  assert (folder/'submission-first.json').exists()
 assert shutil.disk_usage(root).free>15*1024**3
 selected='0' if phase=='first' else '1-2%2'
 job=subprocess.check_output(['sbatch','--parsable','--array='+selected,str(here/'run.sbatch'),receipt],cwd=root,text=True).strip()
 state={'status':'SUBMITTED','phase':phase,'gpu_array_job':job,'tasks':selected,'receipt':receipt,'duration':10,'max_gpu_hours_per_case':1.5,'automatic_extension':False}
 record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

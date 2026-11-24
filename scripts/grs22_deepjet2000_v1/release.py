"""Idempotent release of the single authorized GRS-oriented case after preflight."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys
root=Path(__file__).resolve().parents[2]
here=Path(__file__).resolve().parent
folder=root/'runs/grs22-deepjet2000-v1';folder.mkdir(parents=True,exist_ok=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 record=folder/'submission.json'
 if record.exists():print(record.read_text());sys.exit(0)
 receipt=sys.argv[1]
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 assert shutil.disk_usage(root).free>15*1024**3
 job=subprocess.check_output(['sbatch','--parsable',str(here/'run.sbatch'),receipt],cwd=root,text=True).strip()
 state={'status':'SUBMITTED','gpu_job':job,'receipt':receipt,'fresh_initial_condition':True,'latitude':-22,'target_time':2000,'gpu_release_gate':'exact-driver restart check then full-grid t0-10 pilot; restart unchanged to2000 only if pilot passes','automatic_further_extension':False}
 record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

"""Idempotent release after CPU preflight and completed-donor health check."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
folder=root/'runs/chan-theta-diagnostic4000-v1';folder.mkdir(parents=True,exist_ok=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 record=folder/'submission.json'
 if record.exists():print(record.read_text());sys.exit(0)
 receipt,health=sys.argv[1:3]
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 subprocess.run([str(root/'tools/point4-analysis-env/bin/python'),str(here/'health.py'),'verify',health],cwd=root,check=True)
 assert shutil.disk_usage(root).free>15*1024**3
 job=subprocess.check_output(['sbatch','--parsable',str(here/'run.sbatch'),receipt,health],cwd=root,text=True).strip()
 state={'status':'SUBMITTED','gpu_job':job,'receipt':receipt,'health_receipt':health,'start_time':2000,'target_time':4000,'max_gpu_hours':24,'automatic_extension':False}
 record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

"""Idempotent release of short GPU implementation validation only. No production dispatch."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
folder=root/'runs/chan-total-energy-v1';folder.mkdir(parents=True,exist_ok=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 record=folder/'submission.json'
 if record.exists():print(record.read_text());sys.exit(0)
 receipt=sys.argv[1]
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 assert shutil.disk_usage(root).free>15*1024**3
 job=subprocess.check_output(['sbatch','--parsable',str(here/'run.sbatch'),receipt],cwd=root,text=True).strip()
 state={'status':'SUBMITTED','gpu_validation_job':job,'receipt':receipt,'cases':['compressible theta','compressible total energy'],'grid':[16,16,24],'start_time':0,'target_time':0.04,'omega':3,'beta':0,'jet':False,'max_gpu_minutes':20,'production_submitted':False,'automatic_extension':False}
 record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

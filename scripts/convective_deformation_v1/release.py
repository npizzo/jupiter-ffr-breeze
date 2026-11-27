"""Release only the two user-authorized deformation experiments, once each."""
from pathlib import Path
import fcntl,json,shutil,subprocess,sys
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
folder=root/'runs/convective-deformation-v1';folder.mkdir(parents=True,exist_ok=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 receipt=sys.argv[1]
 subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
 assert shutil.disk_usage(root).free>15*1024**3
 for label,ratio in [('weak',0.5),('strong',2.0)]:
  record=folder/(label+'-submission.json')
  if record.exists():print(record.read_text());continue
  job=subprocess.check_output(['sbatch','--parsable',str(here/'run.sbatch'),receipt,str(ratio),label],cwd=root,text=True).strip()
  state={'status':'SUBMITTED','gpu_job':job,'receipt':receipt,'case':label,'L_over_Ld':ratio,'grid':[256,256,80],'target_turnovers':10,'max_gpu_hours':1,'automatic_extension':False}
  record.write_text(json.dumps(state,indent=2)+'\n');print(json.dumps(state),flush=True)

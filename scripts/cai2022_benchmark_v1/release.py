from pathlib import Path
import fcntl,json,subprocess,sys,shutil,tomllib
root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
folder=root/'runs/cai2022-benchmark-v1';folder.mkdir(parents=True,exist_ok=True)
mode,receipt=sys.argv[1:3]
subprocess.run(['python',str(here/'gate.py'),'verify',receipt],cwd=root,check=True)
with (folder/'launch.lock').open('w') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 if mode=='validation':
  record=folder/'validation-submission.json'
  if record.exists():print(record.read_text());sys.exit(0)
  state={'status':'SUBMITTED_VALIDATION_ONLY','receipt':str(Path(receipt).resolve()),'jobs':{},'production_released':False}
  record.write_text(json.dumps(state,indent=2))
  for case in ['B0','B7']:
   job=subprocess.check_output(['sbatch','--parsable','--job-name=cai-'+case+'-check',str(here/'validate.sbatch'),receipt,case],cwd=root,text=True).strip()
   state['jobs'][case]=job;record.write_text(json.dumps(state,indent=2))
  print(json.dumps(state))
 elif mode=='production':
  case,gpu_path=sys.argv[3:5]
  assert case in ['B0','B7']
  gpu=tomllib.loads(Path(gpu_path).read_text());assert gpu['status']=='PASS'
  if case=='B7':
   cal=json.loads((folder/'B0_CALIBRATION_PASS.json').read_text());assert cal['status']=='PASS'
  record=folder/(case+'-production.json')
  if record.exists():print(record.read_text());sys.exit(0)
  assert shutil.disk_usage('/scratch3/workspace/nicholas_pizzo_uri_edu-jupiter-cai-benchmark-20261127').free>100*1024**3
  job=subprocess.check_output(['sbatch','--parsable','--job-name=cai-'+case+'-production',str(here/'run.sbatch'),receipt,case],cwd=root,text=True).strip()
  state={'status':'SUBMITTED','case':case,'job':job,'receipt':receipt,'gpu_receipt':gpu_path,'automatic_extension':False,'maximum_allocation_hours':24}
  record.write_text(json.dumps(state,indent=2));print(json.dumps(state))

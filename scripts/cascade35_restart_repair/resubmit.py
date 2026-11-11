"""Resume only the failed 35 N chain after a verified diagnostic-check repair."""
from pathlib import Path
import json,subprocess,tomllib,datetime,csv,shutil
r=Path.cwd();record=r/'data/cascade35-restart-v2.json'
assert not record.exists(),'Already submitted; inspect record rather than duplicate jobs'
assert tomllib.loads((r/'data/cascade35-restart-audit-v1/regression.toml').read_text())['status']=='PASS'
old=json.loads((r/'data/cuda122-restart.json').read_text())
ids=[str(j['job_id']) for j in old['jobs'] if j['latitude']==35]
assert not subprocess.check_output(['squeue','-h','-j',','.join(ids),'-o','%i'],text=True).strip(),'Old 35 N jobs still active'
# Keep the existing, independent 55 N chain; the comparison waits for its day 100.
comparison_job=1447348
state55=subprocess.check_output(['sacct','-n','-X','-j',str(comparison_job),'-o','State'],text=True).strip()
assert state55 in ('PENDING','RUNNING','COMPLETED'),state55
archive=r/'data/cascade35-restart-failed-20261111';assert not archive.exists();archive.mkdir()
for name in ['cascade35-v1','cascade35-preflight-v1','cascade35-300-v1']:
 p=r/'data'/name
 if name!='cascade35-preflight-v1':assert not list(p.rglob('frame_*.jld2')),'Do not replace production output'
 shutil.move(str(p),str(archive/name));p.mkdir()
state={'status':'SUBMITTING','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'latitude':35,'regression_job':1459435,'audit_job':1459099,'comparison_day100_job':comparison_job,'archive':str(archive.relative_to(r)),'jobs':[],'model_equations_and_configuration_changed':False}
def save():record.write_text(json.dumps(state,indent=2)+'\n')
def submit(stage,day,script,args,deps=()):
 cmd=['sbatch','--parsable','--kill-on-invalid-dep=yes','--job-name=ffr-beta35-'+stage+str(day)]
 if deps:cmd+=['--dependency=afterok:'+':'.join(map(str,deps))]
 job=int(subprocess.check_output(cmd+[script,*map(str,args)],text=True).strip().split(';')[0]);state['jobs'].append(dict(stage=stage,end_day=day,job_id=job,afterok=list(deps)));save();return job
save()
prev=submit('validation',1,'slurm/cuda122/cascade35-preflight.sbatch',['data/cascade35-preflight-v1'])
for seg in range(1,11):prev=submit('production',seg*10,'slurm/cuda122/cascade35.sbatch',[seg,'data/cascade35-v1','data/cascade35-preflight-v1'],[prev])
submit('rendering',100,'slurm/cascade35-gallery.sbatch',['data/cascade35-v1'],[prev,comparison_job])
for seg in range(11,31):prev=submit('production',seg*10,'slurm/cuda122/cascade35-extend.sbatch',[seg,'data/cascade35-300-v1','data/cascade35-v1','data/cascade35-preflight-v1'],[prev])
submit('rendering',300,'slurm/cascade35-300-gallery.sbatch',['data/cascade35-300-v1','data/cascade35-v1'],[prev])
for name,extended in [('cascade35-v1',False),('cascade35-300-v1',True)]:
 with (r/'data'/name/'submitted-jobs.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['stage','end_day','job','afterok'])
  for j in state['jobs']:
   if (j['end_day']>100)==extended:w.writerow([j['stage'],j['end_day'],j['job_id'],':'.join(map(str,j['afterok']))])
state['status']='SUBMITTED_FULL_VALIDATION_REQUIRED';save();print(json.dumps(state,indent=2))

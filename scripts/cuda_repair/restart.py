"""Unity-only, once-only restart of the two authorized latitude campaigns."""
from pathlib import Path
import subprocess,json,datetime,shutil,csv,hashlib,tomllib
r=Path.cwd(); record=r/'data/cuda122-restart.json'
assert not record.exists(), 'Restart already recorded; inspect instead of duplicating jobs'
probe_path=r/'data/cuda122-probe-v1/PASS.toml'
if probe_path.exists():
 probe=tomllib.loads(probe_path.read_text());assert probe['status']=='PASS'
 for f,h in probe['hashes'].items():assert hashlib.sha256((r/f).read_bytes()).hexdigest()==h,f
else:
 probe_state=subprocess.check_output(['squeue','-h','-j','1447116','-o','%T'],text=True).strip()
 assert probe_state in ('PENDING','RUNNING'), 'Runtime probe is neither passed nor active'
roots=['data/cascadebeta100-v1','data/cascadebeta100-preflight-v1','data/cascade35-v1','data/cascade35-preflight-v1','data/cascade35-300-v1']
oldjobs=[]
for root in roots:
 p=r/root/'submitted-jobs.tsv'
 if p.exists():oldjobs.extend(row['job'] for row in csv.DictReader(p.open(),delimiter='\t'))
if oldjobs:
 active=subprocess.check_output(['squeue','-h','-j',','.join(oldjobs),'-o','%i'],text=True).strip()
 assert not active, 'Previous chain still active: '+active
for root in roots:
 assert not list((r/root).rglob('frame_*.jld2')), 'Scientific output exists; do not reuse this root'
archive=r/'data/cuda122-failed-attempt-20261111';assert not archive.exists();archive.mkdir()
for root in roots:shutil.move(str(r/root),str(archive/Path(root).name));(r/root).mkdir()
state={'status':'SUBMITTING','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'runtime_probe_job':1447116,'runtime':'12.2','archive':str(archive.relative_to(r)),'jobs':[],'science_source_unchanged':True}
def save():record.write_text(json.dumps(state,indent=2)+'\n')
def submit(stage,lat,day,script,args,dependencies=()):
 cmd=['sbatch','--parsable','--kill-on-invalid-dep=yes']
 if dependencies:cmd+=['--dependency=afterok:'+':'.join(map(str,dependencies))]
 cmd+=['--job-name=ffr-beta'+str(lat)+'-'+stage+str(day),script,*map(str,args)]
 job=int(subprocess.check_output(cmd,text=True).strip().split(';')[0]);state['jobs'].append(dict(stage=stage,latitude=lat,end_day=day,job_id=job,afterok=list(dependencies)));save();return job
save();ends={}
for lat,prefix in [(55,'cascadebeta100'),(35,'cascade35')]:
 root='data/'+prefix+'-v1';pre='data/'+prefix+'-preflight-v1'
 prev=submit('validation',lat,1,'slurm/cuda122/'+prefix+'-preflight.sbatch',[pre],[1447116])
 for seg in range(1,11):prev=submit('production',lat,seg*10,'slurm/cuda122/'+prefix+'.sbatch',[seg,root,pre],[prev])
 ends[lat]=prev
 deps=[prev] if lat==55 else [prev,ends[55]]
 submit('rendering',lat,100,'slurm/'+prefix+'-gallery.sbatch',[root],deps)
prev=ends[35]
for seg in range(11,31):prev=submit('production',35,seg*10,'slurm/cuda122/cascade35-extend.sbatch',[seg,'data/cascade35-300-v1','data/cascade35-v1','data/cascade35-preflight-v1'],[prev])
submit('rendering',35,300,'slurm/cascade35-300-gallery.sbatch',['data/cascade35-300-v1','data/cascade35-v1'],[prev])
for root,lat,extension in [('data/cascadebeta100-v1',55,False),('data/cascade35-v1',35,False),('data/cascade35-300-v1',35,True)]:
 with (r/root/'submitted-jobs.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['stage','end_day','job','afterok'])
  for j in state['jobs']:
   if j['latitude']==lat and (j['end_day']>100)==extension:w.writerow([j['stage'],j['end_day'],j['job_id'],':'.join(map(str,j['afterok']))])
state['status']='SUBMITTED_VALIDATION_REQUIRED';save();print(json.dumps(state,indent=2))

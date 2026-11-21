import os,json,subprocess,fcntl,sys
from pathlib import Path
receipt=Path('data/chan-eddies-v3-cpu-2089424/LAUNCH_PASS.json').resolve()
root=Path('data/chan-eddies-v3-screen-2089476').resolve()
gate=['python','scripts/chan_eddies_v3/gate.py','verify',str(receipt)]
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def submit(p,argv):
 if p.exists():
  data=json.loads(p.read_text());assert data['status']=='SUBMITTED';return data['job_id']
 save(p,dict(status='SUBMITTING',argv=argv))
 result=subprocess.run(['sbatch','--parsable']+argv,text=True,capture_output=True,check=True)
 job=result.stdout.strip().split(';')[0];assert job.isdigit(),result.stdout
 save(p,dict(status='SUBMITTED',job_id=job,argv=argv));return job
with (receipt.parent/'CHAIN.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);subprocess.run(gate,check=True)
 review=submit(receipt.parent/'SUBMITTED_review.json',['--dependency=afterok:2089476','--kill-on-invalid-dep=yes','work/chan-eddies-v3-review.sbatch',str(receipt),str(root)])
 jobs={}
 for case in ['case3','case4']:
  jobs[case]=submit(receipt.parent/('SUBMITTED_'+case+'.json'),['--dependency=afterok:'+review,'--kill-on-invalid-dep=yes','slurm/chan-eddies-v3-'+case+'.sbatch',str(receipt),str(root/'PRODUCTION_PASS.toml')])
 state=dict(status='DEPENDENT_CONTINUATIONS_SUBMITTED',screen_job='2089476',review_job=review,production_jobs=jobs,release_condition='Both short GPU cases pass; source fingerprints verified; checkpoint hashes checked at restart.',allocation_minutes=90)
 save(receipt.parent/'CHAIN.json',state);print(json.dumps(state,indent=2))

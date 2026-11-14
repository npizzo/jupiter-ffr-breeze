"""Requeue unchanged simulations only after the recovered analysis gate passes."""
from pathlib import Path
import json,subprocess,datetime
root=Path('data/fplane-point4-v2');results=Path('rendered/fplane-point4-v2')
assert (root/'preflight/PASS.toml').is_file()
assert json.loads((root/'preflight/RECOVERY_PASS.json').read_text())['status']=='PASS'
smoke=json.loads(Path('rendered/fplane-point4-analysis-validation-v3/late_heated/PASS.json').read_text())
assert smoke['status']=='COMPLETE' and smoke['full_decode_passed']
record_path=root/'resubmission-analysis-v3.json'
assert not record_path.exists(),'Already resubmitted'
cases=Path('scripts/fplane_point4_v2/case-order.txt').read_text().splitlines()
for name in cases:
 assert not (root/name).exists(),f'Existing production output {name}'
 assert not (results/name).exists(),f'Existing review output {name}'
r={'status':'SUBMITTING','compute_location':'Unity','root':str(root),'results':str(results),'case_order':cases,'preflight_gpu_job':1725453,'analysis_gate_job':1737737,'superseded_production_array':1725454,'max_concurrent_gpus':2,'source_directory':'scripts/fplane_point4_v2','analysis_environment':'tools/point4-analysis-env','submitted_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
def save():record_path.write_text(json.dumps(r,indent=2)+'\n')
def submit(args):
 p=subprocess.run(['sbatch','--parsable',*args],text=True,capture_output=True,check=True)
 return int(p.stdout.strip().split(';')[0])
save();n=len(cases)
p=submit([f'--array=0-{n-1}%2','slurm/fplane-point4-v2.sbatch','array',str(root)])
r['production_array_job']=p;r['simulation_jobs']={name:f'{p}_{i}' for i,name in enumerate(cases)};save()
a=submit([f'--array=0-{n-1}%2',f'--dependency=aftercorr:{p}','slurm/fplane-point4-review-v3.sbatch',str(root),str(results)])
r['review_array_job']=a;save()
c=submit([f'--dependency=afterok:{a}','slurm/fplane-point4-review-v3.sbatch','compare',str(root),str(results)])
r['comparison_job']=c;r['status']='VALIDATED_SUBMITTED';save();print(json.dumps(r,indent=2))

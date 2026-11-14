"""Submit the authorized bounded campaign, once the preceding resolution audit finishes.
This only schedules Unity jobs; production is gated and limited to two GPUs.
"""
from pathlib import Path
import sys,json,subprocess,os
root=Path('data/fplane-point4-v2');results=Path('rendered/fplane-point4-v2')
resolution=json.loads(Path('rendered/ferrel35-resolution-v1/PASS.json').read_text())
assert resolution['status']=='COMPLETE' and resolution['full_decode_passed']
assert len(sys.argv)==2 and len(sys.argv[1])==40, 'Pass the verified resolution-gallery publication commit'
root.mkdir(exist_ok=False);results.mkdir(exist_ok=False)
record={'status':'SUBMITTING','compute_location':'Unity','root':str(root),'results':str(results),'resolution_gallery_commit':sys.argv[1],
        'case_order':Path(__file__).with_name('case-order.txt').read_text().splitlines(),'max_concurrent_gpus':2,'allocation_hours':16,'estimated_runtime_gate_hours':12,'supersedes':{'preflight_job':1721081,'production_array':1721082,'reason':'L40S ten-day estimate 10.33 hours exceeded the original seven-hour runtime gate; no production ran'},'simulation_jobs':{}}
def save():
 (root/'submission.json').write_text(json.dumps(record,indent=2)+'\n')
def submit(args):
 p=subprocess.run(['sbatch','--parsable',*args],check=True,text=True,capture_output=True)
 return int(p.stdout.strip().split(';')[0])
save()
gate=submit(['--time=01:00:00','slurm/fplane-point4-v2.sbatch','preflight',str(root)]);record['preflight_job']=gate;save()
n=len(record['case_order']);production=submit([f'--array=0-{n-1}%2',f'--dependency=afterok:{gate}','slurm/fplane-point4-v2.sbatch','array',str(root)])
record['production_array_job']=production;record['simulation_jobs']={name:f'{production}_{i}' for i,name in enumerate(record['case_order'])};save()
review=submit([f'--array=0-{n-1}%2',f'--dependency=aftercorr:{production}','slurm/fplane-point4-review-v2.sbatch',str(root),str(results)])
record['review_array_job']=review;save()
compare=submit([f'--dependency=afterok:{review}','slurm/fplane-point4-review-v2.sbatch','compare',str(root),str(results)])
record['comparison_job']=compare;record['status']='GATED_SUBMITTED';save()
print(json.dumps(record,indent=2))

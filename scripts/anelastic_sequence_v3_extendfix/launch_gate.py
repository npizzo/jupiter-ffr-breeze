#!/usr/bin/env python3
"""Fail-closed preflight, submission, and compute-node verification for this campaign."""
import datetime, fcntl, hashlib, json, os, pathlib, subprocess, sys
ROOT = pathlib.Path(__file__).resolve().parents[2]
HERE = pathlib.Path(__file__).resolve().parent
GPU_BATCH = 'slurm/anelastic-sequence-v3-extendfix.sbatch'
CPU_BATCH = 'slurm/anelastic-sequence-v3-extendfix-cpu.sbatch'
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def write(p, obj):
    p=pathlib.Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    temp=p.with_name(p.name+'.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');temp.replace(p)
def manifest(physical_receipt):
    files=list(HERE.glob('*.jl'))+list(HERE.glob('*.py'))
    for name in [GPU_BATCH,CPU_BATCH,'configs/anelastic-sequence-v3/baseline.toml','configs/anelastic-sequence-v3/continuation.toml','scripts/cuda_repair/environment.sh','Project.toml','Manifest.toml','environments/cascade-cuda122/Project.toml']:
        p=ROOT/name
        if not p.is_file(): raise RuntimeError('Required manifest input missing: '+str(p))
        files.append(p)
    optional={}
    for folder in [ROOT,ROOT/'environments/cascade-cuda122']:
        for name in ['Manifest.toml','LocalPreferences.toml']:
            p=folder/name;optional[str(p.relative_to(ROOT))]=digest(p) if p.is_file() else None
        for p in folder.glob('Manifest-v*.toml'): optional[str(p.relative_to(ROOT))]=digest(p)
    jr=pathlib.Path(os.environ['JULIA_BIN']).resolve()
    env={k:os.environ.get(k,'') for k in ['JULIA_LOAD_PATH','JULIA_DEPOT_PATH','JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS']}
    return {'optional_environment_files':optional,'files':{str(p.relative_to(ROOT)):digest(p) for p in sorted(files)},'physical_receipt':str(pathlib.Path(physical_receipt).resolve()),'physical_receipt_sha256':digest(physical_receipt),'julia_binary':str(jr),'julia_binary_sha256':digest(jr),'environment':env}
def verify(receipt,physical_receipt):
    r=json.loads(pathlib.Path(receipt).read_text())
    if r.get('status')!='PASS': raise RuntimeError('No passing exact-driver receipt')
    if r['manifest']!=manifest(physical_receipt): raise RuntimeError('Tested source/configuration/environment changed; rerun exact-driver preflight')
    if digest(r['driver_result'])!=r['driver_result_sha256']: raise RuntimeError('Preflight result changed')
    print('LAUNCH_RECEIPT_VERIFIED',flush=True)
    return r
os.chdir(ROOT)
mode=sys.argv[1]
if mode=='inspect':
    m=manifest(sys.argv[2]);print('MANIFEST_INPUTS_OK '+str(len(m['files']))+' required files',flush=True)
elif mode=='preflight':
    out=pathlib.Path(sys.argv[2]).resolve();physical=sys.argv[3]
    if out.exists(): raise RuntimeError('Refuse to overwrite preflight output')
    before=manifest(physical)
    command=[os.environ['JULIA_BIN'],'--project=.',str(HERE/'pilot.jl'),'--preflight',str(out),physical]
    status=subprocess.run(command).returncode
    if status:
        write(out/'LAUNCH_FAIL.json',{'status':'FAILED','exit_code':status,'time':utc(),'command_scope':'exact-driver CPU preflight'})
        sys.exit(status)
    if before!=manifest(physical): raise RuntimeError('Files changed during preflight')
    result=out/'DRIVER_PASS.toml'
    if not result.is_file(): raise RuntimeError('Driver did not issue its success receipt')
    write(out/'LAUNCH_PASS.json',{'status':'PASS','time':utc(),'manifest':before,'driver_result':str(result),'driver_result_sha256':digest(result),'cpu_job_id':os.environ.get('SLURM_JOB_ID'),'scope':'CPU driver/diagnostics/output/restart; GPU validation still mandatory'})
    print('EXACT_LAUNCH_GATE_PASS '+str(out/'LAUNCH_PASS.json'),flush=True)
elif mode=='verify':
    verify(sys.argv[2],sys.argv[3])
elif mode=='submit':
    receipt=pathlib.Path(sys.argv[2]).resolve();physical=sys.argv[3];verify(receipt,physical)
    state=receipt.with_name('GPU_SUBMISSION.json')
    with open(str(state)+'.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if state.exists():
            previous=json.loads(state.read_text())
            if 'job_id' in previous: print('ALREADY_SUBMITTED '+previous['job_id']);sys.exit(0)
            raise RuntimeError('An earlier submission attempt is unresolved; inspect Slurm before retrying')
        write(state,{'status':'SUBMITTING','time':utc(),'receipt':str(receipt)})
        process=subprocess.run(['sbatch','--parsable',GPU_BATCH,str(receipt),physical],text=True,capture_output=True)
        if process.returncode:
            write(state,{'status':'SUBMISSION_FAILED','time':utc(),'stderr':process.stderr})
            raise RuntimeError(process.stderr)
        job=process.stdout.strip().split(';')[0]
        if not job.isdigit(): raise RuntimeError('Unrecognized sbatch output; inspect before retrying: '+process.stdout)
        write(state,{'status':'SUBMITTED','job_id':job,'time':utc(),'receipt':str(receipt)})
        print('GPU_SUBMITTED '+job,flush=True)
elif mode=='outcome':
    out=pathlib.Path(sys.argv[2]);code=int(sys.argv[3])
    write(out/'PROCESS_OUTCOME.json',{'status':'EXITED_OK_REVIEW_REQUIRED' if code==0 else 'FAILED','exit_code':code,'job_id':os.environ.get('SLURM_JOB_ID'),'time':utc(),'gpu_validation_pass':(out/'VALIDATION_PASS.toml').is_file(),'calibration_review_ready':(out/'REVIEW_REQUIRED.toml').is_file()})
else: raise RuntimeError('Unknown gate mode')

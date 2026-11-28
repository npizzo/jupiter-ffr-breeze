"""One sequential, checkpointed Unity chain to the already authorized t=5000."""
import contextlib, fcntl, hashlib, json, os, subprocess, sys, tempfile, tomllib
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CONTROL=ROOT/'runs/cai-anelastic-continue-v1'
STATE=CONTROL/'chain-state.json'
TARGET=5000.0

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');os.replace(tmp,path)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

@contextlib.contextmanager
def locked():
    CONTROL.mkdir(parents=True,exist_ok=True)
    with (CONTROL/'chain.lock').open('a') as f:
        fcntl.flock(f,fcntl.LOCK_EX);yield

def verify_receipt(receipt):
    subprocess.run([sys.executable,str(HERE/'gate.py'),'verify',str(receipt)],cwd=ROOT,check=True)

def request(checkpoint,time,iteration,number):
    checkpoint=Path(checkpoint).resolve();assert checkpoint.is_file()
    assert 0<=time<TARGET and iteration>0
    data={'checkpoint':str(checkpoint),'restart_time':float(time),'restart_iteration':int(iteration),
          'target_time':TARGET,'checkpoint_sha256':sha(checkpoint),'segment':number}
    p=CONTROL/f'request-{number:04d}.toml';assert not p.exists()
    p.write_text('\n'.join(f'{k} = {json.dumps(v)}' for k,v in data.items())+'\n')
    return p

def submit(state,req,receipt,parent=None):
    # Record before RPC: an ambiguous sbatch timeout must not silently duplicate a job.
    state.update(status='SUBMITTING',request=str(req),receipt=str(receipt))
    write(STATE,state)
    cmd=['sbatch','--parsable']
    if parent:cmd+=['--dependency=afterok:'+str(parent)]
    cmd += [str(HERE/'run.sbatch'),str(receipt),str(req)]
    result=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,timeout=45)
    if result.returncode:
        state.update(status='SUBMISSION_FAILED',error=result.stderr[-2000:]);write(STATE,state)
        raise RuntimeError(result.stderr)
    job=result.stdout.strip().split(';')[0];assert job.isdigit()
    state.setdefault('jobs',[]).append({'job':job,'request':str(req),'dependency':parent})
    state.update(status='QUEUED',current_job=job);write(STATE,state)
    print('CONTINUATION_SUBMITTED',job,'request',req,flush=True)

def outcome_decision(r):
    assert r['target_time']==TARGET
    assert r['status'] in ('COMPLETE','WALL_LIMIT'), 'Numerical failure: no continuation'
    assert r['final_time']>r['start_time'], 'No progress: no continuation'
    assert r['final_iteration']>0 and Path(r['checkpoint']).is_file()
    if r['final_time']>=TARGET-1e-8:
        assert r['status']=='COMPLETE';return 'COMPLETE'
    assert r['status']=='WALL_LIMIT';return 'CONTINUE'

def selftest():
    with tempfile.TemporaryDirectory() as d:
        cp=Path(d)/'checkpoint';cp.write_bytes(b'test')
        r={'target_time':TARGET,'start_time':46.,'final_time':100.,'final_iteration':1000,'checkpoint':str(cp),'status':'WALL_LIMIT'}
        assert outcome_decision(r)=='CONTINUE'
        assert outcome_decision(dict(r,final_time=TARGET,status='COMPLETE'))=='COMPLETE'
        for changes in ({'status':'FAILED'},{'final_time':46.},{'checkpoint':str(cp)+'missing'},{'target_time':6000.}):
            try:outcome_decision(dict(r,**changes))
            except AssertionError:pass
            else:raise AssertionError('Invalid continuation was accepted')
    print('CONTINUATION_CONTROLLER_SELFTEST_PASS',flush=True)

def main():
    mode=sys.argv[1]
    if mode=='selftest':selftest();return
    if mode=='verify-request':
        r=tomllib.loads(Path(sys.argv[2]).read_text());assert r['target_time']==TARGET
        assert sha(r['checkpoint'])==r['checkpoint_sha256']
        print('CHECKPOINT_REQUEST_VERIFIED',r['restart_time'],r['restart_iteration'],flush=True);return
    with locked():
        if mode=='start':
            if STATE.exists():
                print('EXISTING_CHAIN_NO_DUPLICATE',STATE.read_text());return
            receipt=Path(sys.argv[2]).resolve();verify_receipt(receipt)
            req=request(sys.argv[3],float(sys.argv[4]),int(sys.argv[5]),1)
            submit({'target_time':TARGET,'jobs':[],'last_completed_time':float(sys.argv[4])},req,receipt)
        elif mode=='running':
            state=json.loads(STATE.read_text());assert state['current_job']==sys.argv[3]
            state.update(status='STARTING',output_directory=sys.argv[2]);write(STATE,state)
        elif mode=='finish':
            out=Path(sys.argv[2]);receipt=Path(sys.argv[3]);job=os.environ['SLURM_JOB_ID']
            state=json.loads(STATE.read_text());assert state['current_job']==job
            verify_receipt(receipt)
            r=tomllib.loads((out/'OUTCOME.toml').read_text());decision=outcome_decision(r)
            assert not (out/'parent/GUARD_STOP.toml').exists()
            assert Path(r['checkpoint']).resolve().parent==out.resolve()/'parent/checkpoints'
            state['jobs'][-1].update(outcome=r['status'],final_time=r['final_time'],output_directory=str(out))
            state.update(last_completed_time=r['final_time'],checkpoint=r['checkpoint'])
            if decision=='COMPLETE':
                state.update(status='COMPLETE');write(STATE,state);print('TARGET_5000_COMPLETE',flush=True)
            else:
                req=request(r['checkpoint'],r['final_time'],r['final_iteration'],len(state['jobs'])+1)
                submit(state,req,receipt,parent=job)
        elif mode=='failed':
            state=json.loads(STATE.read_text())
            if state['current_job']==sys.argv[3]:
                state.update(status='FAILED_NO_AUTOMATIC_RETRY',output_directory=sys.argv[2],exit_code=sys.argv[4]);write(STATE,state)
        else:raise ValueError(mode)

if __name__=='__main__':main()

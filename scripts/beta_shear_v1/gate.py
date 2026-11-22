import os,sys,json,hashlib,subprocess,fcntl,csv,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent;os.chdir(ROOT)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,x):Path(p).parent.mkdir(parents=True,exist_ok=True);Path(p).write_text(json.dumps(x,indent=2)+'\n')
def manifest():
 files=list(HERE.glob('*.jl'))+list(HERE.glob('*.py'))
 files+=list((ROOT/'configs').glob('beta-shear-v1-*.toml'))+list((ROOT/'slurm').glob('beta-shear-v1-*.sbatch'))
 files += [ROOT/p for p in ['scripts/rotation_channel_v1/model.jl','scripts/rotation_channel_v1/initialization.jl','scripts/cuda_repair/environment.sh','local.env','Project.toml','Manifest.toml','environments/cascade-cuda122/Project.toml']]
 for p in [ROOT,ROOT/'environments/cascade-cuda122']:files+=list(p.glob('Manifest*.toml'))+list(p.glob('LocalPreferences.toml'))
 return {'files':{str(p.relative_to(ROOT)):sha(p) for p in files},'julia':sha(os.environ['JULIA_BIN']),'env':{k:os.environ.get(k,'') for k in ['JULIA_LOAD_PATH','JULIA_DEPOT_PATH','JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS']}}
def verify(receipt):
 r=json.loads(Path(receipt).read_text());assert r['status']=='PASS' and r['manifest']==manifest();assert sha(r['result'])==r['result_sha256'];print('EXACT_DRIVER_RECEIPT_VERIFIED',flush=True)
def submit(receipt,case,stability=None):
 verify(receipt);assert case in ['pilot','shear','control']
 if case!='pilot':
  assert stability;g=json.loads(Path(stability).read_text());assert g['status']=='PASS' and g['manifest']==manifest();assert sha(g['diagnostics'])==g['diagnostics_sha256']
 state=Path(receipt).with_name('SUBMITTED_'+case+'.json')
 with open(str(state)+'.lock','a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if state.exists():raise RuntimeError('Already submitted; inspect state before retry')
  write(state,{'status':'SUBMITTING'})
  args=['sbatch','--parsable','slurm/beta-shear-v1-'+('pilot' if case=='pilot' else 'run')+'.sbatch',receipt,case]
  if stability:args.append(stability)
  res=subprocess.run(args,capture_output=True,text=True,check=True);job=res.stdout.strip().split(';')[0];assert job.isdigit();write(state,{'status':'SUBMITTED','job_id':job});print('SUBMITTED',case,job,flush=True)
mode=sys.argv[1]
if mode=='preflight':
 out=Path(sys.argv[2]).resolve();before=manifest();assert not out.exists()
 subprocess.run([os.environ['JULIA_BIN'],'--startup-file=no','-e','function check(x); if x isa Expr; x.head in (:error,:incomplete) && error(string(x)); foreach(check,x.args); end; end; for f in ARGS; check(Meta.parseall(read(f,String))); end']+[str(p) for p in HERE.glob('*.jl')],check=True)
 subprocess.run([os.environ['JULIA_BIN'],'--project=.',str(HERE/'driver.jl'),'preflight',str(out)],check=True)
 assert before==manifest();p=out/'DRIVER_PASS.toml';assert p.is_file();write(out/'LAUNCH_PASS.json',{'status':'PASS','manifest':before,'result':str(p),'result_sha256':sha(p)})
elif mode=='verify':verify(sys.argv[2])
elif mode=='verify-stability':
 verify(sys.argv[2]);g=json.loads(Path(sys.argv[3]).read_text());assert g['status']=='PASS' and g['manifest']==manifest();assert sha(g['diagnostics'])==g['diagnostics_sha256'];print('STABILITY_RECEIPT_VERIFIED')
elif mode=='submit':submit(sys.argv[2],sys.argv[3],sys.argv[4] if len(sys.argv)>4 else None)
elif mode=='release':
 receipt,out=sys.argv[2:4];verify(receipt);out=Path(out);path=out/'parent/diagnostics.csv';rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(path.open())]
 last=rows[-1];late=[r for r in rows if r['time']>=200]
 tests={'completed_250':last['time']>=250,'finite':all(math.isfinite(v) for r in rows for v in r.values()),'bounded_mach':max(r['max_mach'] for r in rows)<.5,'limited_nonzonal_response':max(r['nonzonal_velocity_rms'] for r in late)<.25*.25,'retained_jet_profile':max(r['jet_profile_rms_error'] for r in late)<.15*.25,'mass_conserved':max(abs(r['mass']/rows[0]['mass']-1) for r in rows)<1e-7,'wall_normal_zero':max(r['wall_w_max'] for r in rows)==0}
 g={'status':'PASS' if all(tests.values()) else 'FAIL','tests':tests,'last':last,'manifest':manifest(),'diagnostics':str(path.resolve()),'diagnostics_sha256':sha(path),'scope':'Unheated 3D jet test with active buoyancy and original dissipation; not proof of asymptotic jet stability or frozen stratification.'}
 stability=out/'STABILITY_PASS.json';write(stability,g);print(json.dumps(g['tests']),flush=True)
 assert g['status']=='PASS','Jet gate failed: production withheld; inspect before changing anything'
 for case in ['shear','control']:submit(receipt,case,str(stability.resolve()))
else:raise RuntimeError(mode)

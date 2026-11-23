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
 files=list(HERE.glob('*.jl'))+list(HERE.glob('*.py'))+list(HERE.glob('*.toml'))+list(HERE.glob('*.sbatch'))
 files+=list((ROOT/'configs').glob('beta-shear-v1-*.toml'))+list((ROOT/'slurm').glob('beta-shear-v1-*.sbatch'))
 files += [ROOT/p for p in ['scripts/cuda_repair/environment.sh','local.env','Project.toml','Manifest.toml','environments/cascade-cuda122/Project.toml']]
 for p in [ROOT,ROOT/'environments/cascade-cuda122']:files+=list(p.glob('Manifest*.toml'))+list(p.glob('LocalPreferences.toml'))
 files += list((ROOT/'depot/packages/Oceananigans/EhPjs/src/TurbulenceClosures/turbulence_closure_implementations/Smagorinskys').glob('*.jl'))
 files += list((ROOT/'depot/packages/Breeze/hw2ti/src').rglob('*.jl'))
 return {'files':{str(p.relative_to(ROOT)):sha(p) for p in files},'julia':sha(os.environ['JULIA_BIN']),'env':{k:os.environ.get(k,'') for k in ['JULIA_LOAD_PATH','JULIA_DEPOT_PATH','JULIA_NUM_THREADS','OPENBLAS_NUM_THREADS']}}
def verify(receipt):
 r=json.loads(Path(receipt).read_text());assert r['status']=='PASS' and r['manifest']==manifest();assert sha(r['result'])==r['result_sha256'];print('EXACT_DRIVER_RECEIPT_VERIFIED',flush=True)
mode=sys.argv[1]
if mode=='preflight':
 out=Path(sys.argv[2]).resolve();before=manifest();assert not out.exists()
 subprocess.run([os.environ['JULIA_BIN'],'--startup-file=no','-e','function check(x); if x isa Expr; x.head in (:error,:incomplete) && error(string(x)); foreach(check,x.args); end; end; for f in ARGS; check(Meta.parseall(read(f,String))); end']+[str(p) for p in HERE.glob('*.jl')],check=True)
 subprocess.run([os.environ['JULIA_BIN'],'--pkgimages=no','--project=.',str(HERE/'driver.jl'),'preflight',str(out)],check=True)
 assert before==manifest();p=out/'DRIVER_PASS.toml';assert p.is_file();write(out/'LAUNCH_PASS.json',{'status':'PASS','manifest':before,'result':str(p),'result_sha256':sha(p)})
elif mode=='verify':verify(sys.argv[2])

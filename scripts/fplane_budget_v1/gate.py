"""Independent native-budget and bounded-output gate; run on Unity."""
from pathlib import Path
import sys,json,tomllib,shutil
import numpy as np
import h5py

def check_case(d):
 files=sorted((d/'native_vorticity').glob('budget_*.jld2'));assert files
 errors=[];last=None
 for p in files:
  with h5py.File(p) as h:
   t=float(h['t'][()]);z=np.asarray(h['Z']).T;mask=np.broadcast_to(np.abs((np.arange(z.shape[1])+.5)/z.shape[1]-.5)[None,:]<=.4,z.shape);res=np.asarray(h['cumulative_residual']).T;parts=np.asarray(h['cumulative_parts']).T;origin=np.asarray(h['origin_Z']).T;n=int(h['checks'][()]);third=np.asarray(h['third_moment_residual']).ravel();thirdparts=np.asarray(h['third_moment_parts']).T
   scale=max(float(np.max(np.abs(z-origin))),float(np.max(np.abs(parts))),1e-30);allow=16384*np.finfo(float).eps*max(float(np.max(np.abs(z))),float(np.max(np.abs(origin))),1e-8)*max(1,n)
   assert np.max(np.abs(res[mask]))<=allow+1e-8*scale,(p,'cumulative closure')
   thirdscale=max(float(np.max(np.abs(thirdparts[1]))),float(np.mean(np.abs(z[mask])**3)),1e-30)
   assert abs(third[1])<=1e-8*thirdscale+16384*np.finfo(float).eps*thirdscale*max(1,n),(p,'third moment closure')
  q=d/'reduced'/f'state_{round(t):012d}.npz'
  if q.exists():
   with np.load(q) as r:
    mismatch=float(np.max(np.abs(z-r['Z'])));assert mismatch<=1e-11*max(float(np.max(np.abs(z))),1e-10)+1e-18,(p,'independent Z mismatch',mismatch)
   errors.append(mismatch)
  last={'day':t/86400,'steps':n,'native_residual_max_interior':float(np.max(np.abs(res[mask]))),'native_residual_relative':float(np.max(np.abs(res[mask])))/max(scale,allow),'third_moment_residual_interior':float(third[1])}
 assert last['steps']>0 and errors
 assert not list((d/'fields3d').glob('state_*.jld2')),'Temporary export did not retire'
 return dict(last,independent_Z_max_error=max(errors))
if __name__=='__main__':
 root=Path(sys.argv[1]);p=root/'preflight';gpu=tomllib.loads((p/'GPU_PASS.toml').read_text());assert gpu['status']=='PASS'
 checks={n:check_case(p/n) for n in ['mature512','mature1024']}
 assert shutil.disk_usage(root).free>60*1024**3,'Insufficient headroom for bounded output'
 (p/'diagnostic_gate.json').write_text(json.dumps({'status':'PASS','checks':checks,'compute_location':'Unity'},indent=2))
 shutil.copyfile(p/'GPU_PASS.toml',p/'PASS.toml');print('NATIVE_AND_RECONSTRUCTION_DIAGNOSTIC_GATE_PASS',checks,flush=True)

"""Independent full-domain and third-moment accounting checks; run on Unity."""
from pathlib import Path
import sys,json,tomllib,shutil,subprocess
import numpy as np
import h5py

def read(p):
 with h5py.File(p) as h:
  keys=['t','origin_time','Z','origin_Z','cumulative_parts','cumulative_residual','third_moment_storage','third_moment_parts','third_moment_residual','central_third_storage','central_third_parts','central_third_residual','conditional_edges','conditional_sum','conditional_count','checks','f','max_relative_step_residual','max_pressure_curl','max_mass_divergence']
  return {k:np.asarray(h[k]).T for k in keys}
def check_case(d):
 files=sorted((d/'native_vorticity').glob('budget_*.jld2'));assert files
 errors=[];last=None
 for p in files:
  h=read(p);z=h['Z'];origin=h['origin_Z'];parts=h['cumulative_parts'];n=int(h['checks']);res=z-origin-parts.sum(2)
  assert np.allclose(res,h['cumulative_residual'],rtol=0,atol=1e-20)
  scale=max(float(abs(z-origin).max()),float(abs(parts).max()),1e-30)
  allowance=16384*np.finfo(float).eps*max(float(abs(z).max()),float(abs(origin).max()),1e-8)*max(1,n)
  assert abs(res).max()<=allowance+1e-8*scale,(p,'full closure',abs(res).max(),scale)
  interior=np.broadcast_to(abs((np.arange(z.shape[1])+.5)/z.shape[1]-.5)[None,:]<=.4,z.shape)
  for region,mask in enumerate([np.ones(z.shape,bool),interior,~interior]):
   for prefix in ['third_moment','central_third']:
    q=z[mask];o=origin[mask]
    if prefix=='central_third':q=q-q.mean();o=o-o.mean()
    actual=np.mean(q**3)-np.mean(o**3);terms=h[prefix+'_parts'][region];reported=float(h[prefix+'_storage'][region])
    mscale=max(float(np.mean(abs(q)**3)),float(np.max(abs(terms))),1e-40)
    assert abs(actual-reported)<1e-10*mscale+1e-40
    assert abs(actual-terms.sum())<1e-8*mscale+16384*np.finfo(float).eps*mscale*max(1,n),(p,prefix,region)
  q=d/'reduced'/f'state_{round(float(h["t"])):012d}.npz'
  if q.exists():
   with np.load(q) as r:
    mismatch=float(abs(z-r['Z']).max());assert mismatch<1e-11*max(float(abs(z).max()),1e-10)+1e-18
   errors.append(mismatch)
  last={'day':float(h['t'])/86400,'steps':n,'max_step_relative_residual_full':float(h['max_relative_step_residual']),'cumulative_residual_max_full':float(abs(res).max()),'pressure_curl_max':float(h['max_pressure_curl']),'depth_mean_divergence_max':float(h['max_mass_divergence'])}
 assert last['steps']>0 and errors
 last['independent_Z_max_error']=max(errors)
 return last
if __name__=='__main__':
 root=Path(sys.argv[1]);p=root/'preflight';gpu=tomllib.loads((p/'GPU_PASS.toml').read_text());assert gpu['status']=='PASS'
 # Finished preflight fields are temporary validation exports; reduce after exact paired checks.
 for name in ['baseline','plus','mirror']:
  for field in sorted((p/name/'fields3d').glob('state_*.jld2')):
   subprocess.run([sys.executable,str(Path(__file__).with_name('reduce.py')),str(field),str(p/name/'reduced')],check=True)
 checks={name:check_case(p/name) for name in ['plus','mirror']}
 a=read(sorted((p/'plus'/'native_vorticity').glob('budget_*.jld2'))[-1]);b=read(sorted((p/'mirror'/'native_vorticity').glob('budget_*.jld2'))[-1])
 assert np.max(abs(b['Z']+a['Z'][::-1]))<1e-6*max(np.max(abs(a['Z'])),1e-10)
 assert np.max(abs(b['cumulative_parts']+a['cumulative_parts'][::-1]))<1e-5*max(np.max(abs(a['cumulative_parts'])),1e-15)
 assert shutil.disk_usage(root).free>30*1024**3
 (p/'diagnostic_gate.json').write_text(json.dumps({'status':'PASS','checks':checks,'compute_location':'Unity'},indent=2))
 shutil.copyfile(p/'GPU_PASS.toml',p/'PASS.toml')
 print('SIGN_FULL_DOMAIN_GATE_PASS',checks,flush=True)

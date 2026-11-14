"""Analytic and saved-field checks executed on Unity before production release."""
from pathlib import Path
import sys,json,tempfile,shutil
import numpy as np,h5py
from diagnostics import poisson,lap,field_analysis,fit_relation,moments,smooth
root=Path(sys.argv[1]);checks=[]
assert shutil.disk_usage(root).free>80*1024**3, 'Less than 80 GiB remains for the bounded full-field campaign'
nx,ny,nz=32,24,8;hx,hy,hz=20000.,25000.,1500.
x=(np.arange(nx)+.5)*hx;y=(np.arange(ny)+.5)*hy
p=np.sin(2*np.pi*x[:,None]/(nx*hx))*np.sin(2*np.pi*y[None,:]/(ny*hy))
z=lap(p,hx,hy);sol=poisson(z,hx,hy)
assert np.max(abs(sol-p))<1e-12;checks.append('Periodic-x, channel-y Poisson inversion and sign')
assert moments(np.array([-1.,0.,1.]))['raw_skewness']==0;checks.append('Symmetric vorticity has zero skewness')
assert np.allclose(smooth(np.zeros((nx,ny)),100,hx,hy),0);checks.append('Wall-aware filter preserves zero')
q=np.linspace(-2,2,1000);fit=fit_relation(q,1+2*q+.4*q*q)
assert fit['quadratic_R2']>1-1e-12 and abs(fit['raw_c_over_b']-.2)<1e-10;checks.append('Quadratic dimensional coefficient conversion')
# Native-face velocity from a corner streamfunction; density changes with height.
xf=np.arange(nx)*hx;yf=np.arange(ny+1)*hy
psi=1e6*np.sin(2*np.pi*xf[:,None]/(nx*hx))*np.sin(2*np.pi*yf[None,:]/(ny*hy))
uf=np.broadcast_to((-np.diff(psi,axis=1)/hy)[:,:,None],(nx,ny,nz)).copy()
vf=np.broadcast_to(((np.roll(psi,-1,axis=0)-psi)/hx)[:,:,None],(nx,ny+1,nz)).copy()
wf=np.zeros((nx,ny,nz+1));rho=np.exp(-np.arange(nz)/nz)
path=root/'manufactured_barotropic.h5'
with h5py.File(path,'w') as d:
 for key,a in [('u_face',uf),('v_face',vf),('w_face',wf)]:d[key]=a.T
 for key,a in [('density',rho),('dx',hx),('dy',hy),('dz',hz),('f',2.8e-4),('t',0.)]:d[key]=a
row,m=field_analysis(path)
assert row['native_depth_mean_divergence_max']<1e-16
assert np.max(abs(m['vertical_source']))<1e-25 and np.max(abs(m['horizontal_source']))<1e-25
checks.append('Density-weighted barotropic flow has no baroclinic V or H source')
path.unlink()
# Real GPU snapshots validate orientation, density, native-face continuity, and finiteness.
for case in ['late_heated','late_control','late_restart_reversed','rest_q540','rest_q162_flipf','rest_q162_fine']:
 files=sorted((root/case/'fields3d').glob('state_*.jld2'));assert len(files)==2
 row,m=field_analysis(files[-1])
 assert row['native_depth_mean_divergence_max']<1e-10
 assert all(np.isfinite(m[k]).all() for k in ['Z','Psi','vertical_source','horizontal_source'])
 checks.append(case+': native 3D checkpoint export and finite diagnostics')
(root/'DIAGNOSTICS_PASS.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','checks':checks},indent=2)+'\n')
print('POINT4_ALL_PREFLIGHT_CHECKS_PASS',len(checks),flush=True)

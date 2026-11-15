"""Independent field reduction only; Unity. Preserve scratch on any failure."""
from pathlib import Path
import sys,json
import numpy as np
import h5py
p=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
assert p.parent.name=='fields3d' and p.name.startswith('state_')
def dx(a,h):return (np.roll(a,-1,0)-np.roll(a,1,0))/(2*h)
def dy(a,h,parity=1):
 q=np.concatenate((parity*a[:,0:1],a,parity*a[:,-1:]),axis=1)
 return (q[:,2:]-q[:,:-2])/(2*h)
with h5py.File(p) as d:
 uf,vf,wf=[np.asarray(d[k]).T for k in ['u_face','v_face','w_face']]
 rho=d['density'][:];weights=rho[None,None,:]/rho.sum();h=float(d['dx'][()]);hy=float(d['dy'][()]);t=float(d['t'][()]);f=float(d['f'][()])
 u=(uf+np.roll(uf,-1,0))/2;v=(vf[:,:-1]+vf[:,1:])/2
 zeta=dx(v,h)-dy(u,hy);Z=np.sum(weights*zeta,axis=2)
 U=np.sum(weights*u,axis=2);V=np.sum(weights*v,axis=2)
 up=u-U[:,:,None];vp=v-V[:,:,None]
 rxx=np.sum(weights*up*up,axis=2);rxy=np.sum(weights*up*vp,axis=2);ryy=np.sum(weights*vp*vp,axis=2)
 C=dx(-dx(rxy,h)-dy(ryy,hy),h)-dy(-dx(rxx,h)-dy(rxy,hy,-1),hy)
 assert np.isfinite(Z).all() and np.isfinite(C).all()
np.savez_compressed(out/(p.stem+'.npz'),Z=Z,C=C,U=U,V=V,t=t,f=f)
(out/(p.stem+'.json')).write_text(json.dumps({'time_days':t/86400,'z_rms':float(np.sqrt(np.mean(Z**2))),'shape':list(Z.shape),'compute_location':'Unity'}))
with np.load(out/(p.stem+'.npz')) as d:assert d['Z'].shape==Z.shape and np.isfinite(d['Z']).all()
p.unlink()
print('REDUCED_TEMPORARY_FIELD',p.name,flush=True)

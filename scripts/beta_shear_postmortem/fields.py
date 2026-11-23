"""Read-only analysis, executed on Unity. Arrays are z,y,x in h5py."""
import json, pathlib, h5py, numpy as np
root=pathlib.Path('data/beta-shear-v1-pilot-2134379/parent')
out=pathlib.Path('data/beta-shear-v1-postmortem');out.mkdir(exist_ok=True)
rows=[]
for path in sorted((root/'fields3d').glob('state*.jld2')):
 with h5py.File(path) as h:
  rho=h['density'][:].astype(float);T=h['temperature'][:].astype(float)
  uf=h['u_face'][:].astype(float);vf=h['v_face'][:].astype(float);wf=h['w_face'][:].astype(float)
  t=float(h['t'][()]);z=h['z'][:]
 u=.5*(uf+np.roll(uf,-1,axis=2));v=.5*(vf[:,:-1]+vf[:,1:]);w=.5*(wf[:-1]+wf[1:])
 nz,ny,nx=rho.shape;dy=8/ny;dz=1/nz;dx=4/nx;y=(np.arange(ny)+.5)*dy-4
 mach=np.sqrt(u*u+v*v+w*w)/np.sqrt(5*T/3)
 lnth=.6*np.log(T)-.4*np.log(rho)
 # Ertel PV sign proxy omega_absolute dot grad(log(theta)); rho factor positive.
 def ddx(a):return (np.roll(a,-1,axis=2)-np.roll(a,1,axis=2))/(2*dx)
 def ddy(a):return np.gradient(a,dy,axis=1)
 def ddz(a):return np.gradient(a,dz,axis=0)
 fh=6*np.cos(np.deg2rad(67.5));f=6*np.sin(np.deg2rad(67.5))+.003211688852166031*y[None,:,None]
 pv=(ddy(w)-ddz(v))*ddx(lnth)+(ddz(u)-ddx(w)+fh)*ddy(lnth)+(ddx(v)-ddy(u)+f)*ddz(lnth)
 n2=12*ddz(lnth)
 idx=np.unravel_index(np.argmin(T),T.shape);im=np.unravel_index(np.argmax(mach),T.shape)
 r=dict(time=t,Tmin=float(T[idx]),Tmin_xyz=[float((idx[2]+.5)*dx),float(y[idx[1]]),float(z[idx[0]])],max_mach=float(mach[im]),mach_xyz=[float((im[2]+.5)*dx),float(y[im[1]]),float(z[im[0]])],bands={})
 for name,mask in [('deep',z<.65),('transition',(z>=.65)&(z<.9)),('cap',z>=.9)]:
  zonal=[a.mean(axis=2,keepdims=True) for a in [u,v,w]]
  b=dict(Tmin=float(T[mask].min()),Tmean=float(T[mask].mean()),max_mach=float(mach[mask].max()),u_rms=float(np.sqrt(np.mean(u[mask]**2))),v_rms=float(np.sqrt(np.mean(v[mask]**2))),w_rms=float(np.sqrt(np.mean(w[mask]**2))),v_zonal_rms=float(np.sqrt(np.mean(zonal[1][mask]**2))),w_zonal_rms=float(np.sqrt(np.mean(zonal[2][mask]**2))),eddy_rms=float(np.sqrt(sum(np.mean((a[mask]-az[mask])**2) for a,az in zip([u,v,w],zonal)))),pv_negative_fraction=float(np.mean(pv[mask]<0)),pv_percentiles=np.percentile(pv[mask],[0,10,50,90,100]).tolist(),N2_percentiles=np.percentile(n2[mask],[0,10,50,90,100]).tolist())
  r['bands'][name]=b
 rows.append(r)
 print(json.dumps(r),flush=True)
(out/'field-summary.json').write_text(json.dumps(rows,indent=2)+'\n')

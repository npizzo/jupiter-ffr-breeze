"""Read saved fields on Unity only. Diagnostic Rossby is NOT Cai's reference Ro_c."""
import csv,json,sys,tomllib
from pathlib import Path
import h5py,numpy as np
base=Path.cwd(); out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
cases={'none':'1990851','cyclonic':'1990854','anticyclonic':'1990856'}
rows=[]
for case,job in cases.items():
 root=base/f'data/strong256-direct-v1-{job}'
 cfg=tomllib.loads((root/'day1_50/metadata.toml').read_text())
 f=cfg['f']; Hc=cfg['cap_height']; dx=cfg['lx']/cfg['nx'];dy=cfg['ly']/cfg['ny'];nu=cfg['viscosity']
 files={p.name:p for seg in ('day0_1','day1_50') for p in (root/seg/'fields3d').glob('frame_*.jld2')}
 for name,p in sorted(files.items()):
  with h5py.File(p) as d:
   t=float(d['t'][()]);z=d['z'][:];rho=d['density'][:];uf=d['u_face'][:].astype(float);vf=d['v_face'][:].astype(float);wf=d['w_face'][:].astype(float)
  u=.5*(uf+np.roll(uf,-1,axis=2));v=.5*(vf[:,:-1]+vf[:,1:]);w=.5*(wf[:-1]+wf[1:])
  m=(z>=10000)&(z<=70000);core=slice(int(np.ceil(.1*cfg['ny']))-1,int(np.floor(.9*cfg['ny'])))
  q=rho[m]/rho[m].sum();Uc=(u[m,:, :]-u[m].mean(axis=(1,2),keepdims=True));Vc=v[m]-v[m].mean(axis=(1,2),keepdims=True)
  # Remove depth-coherent horizontal motion to avoid calling an LSV 'convection'.
  Uc-=np.sum(q[:,None,None]*Uc,axis=0,keepdims=True);Vc-=np.sum(q[:,None,None]*Vc,axis=0,keepdims=True)
  W=w[m]-w[m].mean(axis=(1,2),keepdims=True)
  rms=lambda a:float(np.sqrt(np.sum(q*np.mean(a[:,core,:]**2,axis=(1,2)))))
  wrms=rms(W);urms=np.sqrt(rms(Uc)**2+rms(Vc)**2+wrms**2)
  # Channel w is even across walls. Reflection supplies its Neumann spectral extension.
  we=np.concatenate((W,W[:,::-1,:]),axis=1);A=np.fft.rfft2(we,axes=(1,2),norm='forward');E=abs(A)**2
  factor=np.ones(E.shape[-1]);factor[1:-1]=2;E*=factor[None,None,:]
  ky=2*np.pi*np.fft.fftfreq(we.shape[1],dy)[:,None];kx=2*np.pi*np.fft.rfftfreq(we.shape[2],dx)[None,:];K=np.hypot(kx,ky)
  Ek=np.sum(q[:,None,None]*E,axis=0);Ek[0,0]=0;kbar=float(np.sum(K*Ek)/max(Ek.sum(),1e-300));ell=1/max(kbar,1e-300)
  meanu=np.mean(u[m],axis=2);zshear=-np.gradient(meanu,dy,axis=1);central=zshear[:,core]/f
  theta_path=root/'day1_50'/name
  if not theta_path.exists():theta_path=root/'day0_1'/name
  with h5py.File(theta_path) as d:
   hp=d['heat_flux_profile'][:];dp=d['diffusive_heat_flux_profile'][:]
  r=dict(case=case,day=t/86400,w_rms=wrms,baroclinic_velocity_rms=urms,convective_gradient_length_m=ell,convective_wavelength_m=2*np.pi*ell,Ro_w_local=wrms/(abs(f)*ell),Ro_baroclinic_local=urms/(abs(f)*ell),Ro_w_depth=wrms/(abs(f)*Hc),Re_w_depth=wrms*Hc/nu,Re_w_local=wrms*ell/nu,E_vertical_f=nu/(abs(f)*Hc**2),shear_over_f_min=float(central.min()),shear_over_f_max=float(central.max()),shear_over_f_center=float(np.sum(q*zshear[:,cfg['ny']//2])/f),w_energy_shorter_than_6dx=float(Ek[K>2*np.pi/(6*dx)].sum()/max(Ek.sum(),1e-300)),upward_convective_flux_mid=float(hp[len(hp)//2]),upward_diffusive_flux_mid=float(dp[len(dp)//2]))
  rows.append(r)
 print('AUDITED',case,len(files),flush=True)
with (out/'audit.csv').open('w') as h:
 writer=csv.DictWriter(h,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
summary={}
for case in cases:
 rr=[r for r in rows if r['case']==case and 40<=r['day']<=50]
 summary[case]={k:dict(mean=float(np.mean([r[k] for r in rr])),min=float(min(r[k] for r in rr)),max=float(max(r[k] for r in rr))) for k in rr[0] if k not in ('case','day')}
(out/'summary.json').write_text(json.dumps(dict(status='COMPLETE',late_window_days=[40,50],late_statistics=summary,definitions={'ell':'1 / energy-weighted mean horizontal w wavenumber; wavelength=2pi ell','velocity':'depth-varying horizontal residual after removing horizontal means and convective-layer density-weighted mean, plus w; central ~80% of channel','Ro_c_Cai':'UNAVAILABLE: requires a matched nonrotating reference; do not substitute initial random Ro or diagnostic Ro','source':'existing 256x256x64 day0-50 snapshots; no rerun','limitations':['daily snapshots','centered fields','spectrum includes full reflected channel; RMS excludes walls','short-wave metric is resolution warning, not a truncation error estimate','not an inverse-transfer or vorticity-budget calculation']},rows=len(rows)),indent=2)+'\n')
print('REGIME_AUDIT_COMPLETE',flush=True)

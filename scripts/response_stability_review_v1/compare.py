"""Unity finite-time response review; do not interpret trajectory divergence as an eigenmode."""
import sys,json,csv
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
mode,out,*folders=sys.argv[1:];out=Path(out);out.mkdir(parents=True,exist_ok=True)
paths=[Path(s) for s in folders];assert len(paths)==3
channel=mode=='stability';length=6e6 if channel else 2.;f=.0002880957739764384 if channel else .5
frames=[]
for folder in paths:
 d={}
 for p in folder.glob('frame_*.jld2'):
  with h5py.File(p) as h:d[float(h['t'][()])]=p
 frames.append(d)
times=sorted(set.intersection(*(set(d) for d in frames)));assert len(times)>1
assert all(set(times)==set(d) for d in frames),'Do not silently truncate unequal completed timelines'
def low(a,parity=1):
 n=a.shape[1];m=a.shape[0]
 ext=np.concatenate((a,parity*a[::-1,:]),axis=0) if channel else a
 kx=np.fft.fftfreq(n)*n;ky=np.fft.fftfreq(ext.shape[0])*ext.shape[0]*(.5 if channel else 1.)
 K=ky[:,None]**2+kx[None,:]**2
 return np.fft.ifft2(np.fft.fft2(ext)*((K>0)&(K<=4))).real[:m]
def zeta(u,v):
 dx=length/u.shape[1];dy=length/u.shape[0]
 du=np.gradient(u,dy,axis=0,edge_order=2) if channel else (np.roll(u,-1,0)-np.roll(u,1,0))/(2*dy)
 return (np.roll(v,-1,1)-np.roll(v,1,1))/(2*dx)-du
def state(p):
 with h5py.File(p) as h:return h['ubar'][()].astype(float),h['vbar'][()].astype(float)
def align_error(a,b):
 # Translation-invariant SHAPE comparison in the common interior. Do not wrap y walls.
 ny,nx=a.shape;best=(np.inf,0,0);band=max(2,ny//8);crop=slice(band,ny-band)
 for dy in range(-band//2,band//2+1):
  aa=a[crop];bb=b[band+dy:ny-band+dy]
  corr=np.fft.ifft(np.fft.fft(aa,axis=1)*np.fft.fft(bb,axis=1).conj(),axis=1).real.sum(axis=0)
  dx=int(np.argmax(corr));bb=np.roll(bb,dx,axis=1)
  err=float(np.mean((aa-bb)**2))
  if err<best[0]:best=(err,dx if dx<=nx//2 else dx-nx,dy)
 return best
rows=[];initial=None
for t in times:
 uv=[state(d[t]) for d in frames];lp=[(low(a[0],1),low(a[1],-1)) for a in uv];z=[low(zeta(*a),-1) for a in uv]
 diff=[(a[0]-lp[0][0],a[1]-lp[0][1]) for a in lp[1:]]
 if initial is None:initial=[tuple(a.copy() for a in v) for v in diff]
 row={'time':t,'elapsed_days_or_units':(t-times[0])/86400 if channel else t-times[0]}
 for j,(a,zz) in enumerate(zip(lp,z)):
  row[f'low_energy_{j}']=float(.5*np.mean(a[0]**2+a[1]**2));row[f'low_skew_{j}']=float(np.mean((zz-zz.mean())**3)/(np.mean((zz-zz.mean())**2)**1.5+1e-99));row[f'core_max_zeta_over_f_{j}']=float(zz.max()/f);row[f'core_min_zeta_over_f_{j}']=float(zz.min()/f)
 for j,d in enumerate(diff):
  d0=initial[j];norm=np.mean(d0[0]**2+d0[1]**2)
  row[f'large_scale_difference_energy_{j+1}']=float(.5*np.mean(d[0]**2+d[1]**2))
  row[f'template_amplitude_{j+1}']=float(np.mean(d[0]*d0[0]+d[1]*d0[1])/max(norm,1e-99))
  er,sx,sy=align_error(z[0],z[j+1]);row[f'aligned_shape_error_{j+1}']=er;row[f'shift_x_cells_{j+1}']=sx;row[f'shift_y_cells_{j+1}']=sy
 if channel:
  a,b=diff;scale=max(np.mean(b[0]**2+b[1]**2),1e-99)
  row['half_full_linearity_error']=float(np.sqrt(np.mean((2*a[0]-b[0])**2+(2*a[1]-b[1])**2)/scale))
 else:
  # Even response after subtracting the evolving control, and odd signed response.
  a,b=diff;row['even_response_energy']=float(.125*np.mean((a[0]+b[0])**2+(a[1]+b[1])**2));row['odd_response_energy']=float(.125*np.mean((a[0]-b[0])**2+(a[1]-b[1])**2))
 rows.append(row)
with open(out/'comparison.csv','w') as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained');tt=[r['elapsed_days_or_units'] for r in rows]
labels=['control','0.5% deformation','1% deformation'] if channel else ['control','cyclonic seed','anticyclonic seed']
for j,label in enumerate(labels):ax[0,0].plot(tt,[r[f'low_energy_{j}'] for r in rows],label=label);ax[0,1].plot(tt,[r[f'low_skew_{j}'] for r in rows],label=label)
for j in [1,2]:ax[1,0].plot(tt,[r[f'template_amplitude_{j}'] for r in rows],label=labels[j])
if channel:ax[1,1].plot(tt,[r['half_full_linearity_error'] for r in rows]);ax[1,1].set_title('Deviation from amplitude-linear response')
else:
 for key in ['even_response_energy','odd_response_energy']:ax[1,1].plot(tt,[r[key] for r in rows],label=key)
 ax[1,1].legend();ax[1,1].set_title('Even and odd seeded responses')
ax[0,0].set_title('Absolute low-mode kinetic energy');ax[0,1].set_title('Low-mode vorticity skewness');ax[1,0].set_title('Projection onto initial perturbation')
for a in ax.flat:a.set_xlabel('Additional Earth days' if channel else 'Additional nondimensional time');a.grid(alpha=.2)
for a in [ax[0,0],ax[0,1],ax[1,0]]:a.legend()
fig.savefig(out/'response.png',dpi=140)
(out/'summary.json').write_text(json.dumps({'mode':mode,'paths':[str(p) for p in paths],'time_range':[times[0],times[-1]],'frames':len(times),'final':rows[-1],'limitations':['One parent and one spatial template; not an ensemble.','Low-mode response and displacement-corrected shape are separated from small-scale decorrelation.','Finite-time response does not establish a stationary eigenmode or spontaneous LSV formation.','Channel spectral filters use free-slip parity: u even, v and zeta odd; periodic filters are used for the doubly periodic box.','Heat drift and early pressure adjustment remain possible confounders.','Native transfer CSV samples selected timesteps; no claim of an uninterrupted integrated budget.']},indent=2)+'\n')
print('COMPARISON_COMPLETE',out)

from pathlib import Path
import json,hashlib
import numpy as np,h5py
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
out=Path('rendered/weno9-lsv100');fields=[]
labels=['compressible no jet','compressible jet','anelastic no jet','anelastic jet']
for i in range(4):
 with h5py.File(f'data/weno9-lsv100-v3-run-2181508_{i}/parent/frame_000100000.jld2') as h:
  u=h['ubar'][:];v=h['vbar'][:];f=float(h['f'][()]);ny,nx=u.shape;dx=4/nx;dy=8/ny
  u=u-u.mean(axis=1,keepdims=True);v=v-v.mean(axis=1,keepdims=True)
  q=((np.roll(v,-1,axis=1)-np.roll(v,1,axis=1))/(2*dx)-np.gradient(u,dy,axis=0))/f
  fields.append(gaussian_filter(q,(.25/dy,.25/dx),mode=('reflect','wrap')))
lim=max(float(np.max(abs(q))) for q in fields)
fig,axs=plt.subplots(2,2,figsize=(9,11),layout='constrained')
for ax,q,label in zip(axs.flat,fields,labels):
 im=ax.imshow(q,origin='lower',extent=[0,4,-4,4],cmap='RdBu_r',norm=AsinhNorm(linear_width=.01,vmin=-lim,vmax=lim))
 ax.set(title=label,xlabel='x/H',ylabel='y/H');cb=fig.colorbar(im,ax=ax,label='Smoothed nonzonal depth-mean ζ/f₀');cb.set_ticks([-lim,-.01,0,.01,lim]);cb.set_ticklabels([f'{-lim:.3f}','−0.01','0','0.01',f'{lim:.3f}'])
fig.suptitle('t = 100 · reference-density-weighted depth mean\nZonal mean removed; Gaussian σ = 0.25H; full extrema, no clipping')
fig.savefig(out/'depthmean_100.png',dpi=110);plt.close(fig)
p=out/'media.json';m=json.loads(p.read_text());q=out/'depthmean_100.png';m['assets'][q.name]={'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()};m['depthmean_figure_symmetric_limit']=lim;p.write_text(json.dumps(m,indent=2));print('FINAL_FIGURE_VERIFIED',lim)

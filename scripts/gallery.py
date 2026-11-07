"""Render actual saved solver output. Never fabricates demonstration results."""
import argparse, json, pathlib, tomllib
import numpy as np, h5py, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
p=argparse.ArgumentParser();p.add_argument('run');p.add_argument('--output',default='gallery');a=p.parse_args()
r=pathlib.Path(a.run); dest=pathlib.Path(a.output)/r.name;dest.mkdir(parents=True,exist_ok=True)
m=tomllib.loads((r/'metadata.toml').read_text());files=sorted(r.glob('frame_*.jld2'))
if not files: raise SystemExit('No solver frames found')
frames=[]
for file in files:
 with h5py.File(file) as d:
  theta=d['theta'][:].T;w=d['w'][:].T;vort=d['zeta'][:].T;wxz=d['w_xz'][:].T;t=float(d['t'][()])
 fig,axes=plt.subplots(2,2,figsize=(10,8),layout='constrained')
 fields=[(theta-theta.mean(),'Potential-temperature anomaly [K]'),(vort/m['f'],'Vertical vorticity / f'),(w,'Vertical velocity [m/s]'),(wxz,'Vertical velocity: vertical section [m/s]')]
 for j,(ax,(v,title)) in enumerate(zip(axes.flat,fields)):
  lim=max(float(np.max(np.abs(v))),1e-12)
  extent=[0,m['lx']/1e3,0,(m['height'] if j==3 else m['ly'])/1e3]
  im=ax.imshow(v.T,origin='lower',extent=extent,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim)
  ax.set(title=title,xlabel='x [km]',ylabel='z [km]' if j==3 else 'y [km]');fig.colorbar(im,ax=ax,shrink=.8)
 fig.suptitle(f"{m['name']} | t = {t/3600:.2f} h | dry anelastic pilot\n{m['nx']} × {m['ny']} × {m['nz']} | {m['architecture']} | autoscaled panels")
 fn=dest/(file.stem+'.png');fig.savefig(fn,dpi=110);plt.close(fig);frames.append(fn)
Image.open(frames[-1]).convert('RGB').resize((750,600),Image.Resampling.LANCZOS).quantize(colors=64,dither=Image.Dither.NONE).save(dest/'latest.png',optimize=True)
# A compact animation preview; retain the full-resolution final PNG and frames.
stride=4 if len(frames)>12 else 1
preview_frames=frames[::stride]
if preview_frames[-1] != frames[-1]: preview_frames.append(frames[-1])
ims=[Image.open(f).convert('RGB').resize((400,320),Image.Resampling.LANCZOS).quantize(colors=32,dither=Image.Dither.NONE) for f in preview_frames]
ims[0].save(dest/'evolution.gif',save_all=True,append_images=ims[1:],duration=200*stride,loop=0)
(dest/'metadata.json').write_text(json.dumps(m,indent=2)+'\n')
(dest/'README.md').write_text(f"# {m['name']}\n\nActual solver output. Status: **{m['status']}**. This is an engineering / convection test, not a demonstrated FFR.\n\n![Final snapshot](latest.png)\n\n![Evolution](evolution.gif)\n\nPanels autoscale independently in each frame; compare amplitudes using the color bars. [Run parameters](metadata.json).\n")
print(dest)

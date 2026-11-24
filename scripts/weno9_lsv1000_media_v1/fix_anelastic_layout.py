from pathlib import Path
import json,hashlib,subprocess
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter
out=Path('rendered/weno9-lsv1000');m=json.loads((out/'media.json').read_text());limits=m['vorticity_limits'];heights=m['heights'];frames={}
for i,label in [(3,'jet'),(2,'nojet')]:
 frames[label]={}
 for d in [Path(f'data/weno9-startup-v2-pilot-2177387_{i}/parent'),Path(f'data/weno9-lsv100-v3-run-2181508_{i}/parent'),Path(f'data/weno9-lsv1000-v4-run-2187099_{i}/parent')]:
  for p in d.glob('frame_*.jld2'):
   t=int(p.stem.removeprefix('frame_'))/1000
   if abs(t-round(t))<1e-7 and round(t)%5==0:frames[label][round(t)]=p
 assert sorted(frames[label])==list(range(0,1001,5))
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':9})
fig,axs=plt.subplots(4,3,figsize=(14,15),dpi=90)
fig.subplots_adjust(left=.06,right=.935,bottom=.045,top=.925,wspace=.60,hspace=.52)
images=[];normtr=FuncNorm((lambda x:np.arcsinh((x-.1)/.005),lambda x:.1+.005*np.sinh(x)),vmin=0,vmax=1)
for row in range(4):
 for j in range(3):
  tracer=row%2==1;norm=normtr if tracer else AsinhNorm(linear_width=.01,vmin=-limits[j],vmax=limits[j])
  im=axs[row,j].imshow(np.zeros((256,128)),origin='lower',extent=[0,4,-4,4],aspect='auto',interpolation='nearest',cmap='BrBG' if tracer else 'RdBu_r',norm=norm);images.append(im)
  cb=fig.colorbar(im,ax=axs[row,j],fraction=.045,pad=.035)
  if tracer:cb.set_ticks([0,.05,.09,.1,.11,.15,1],labels=['0','.05','.09','.10','.11','.15','1']);cb.set_label('Tracer C')
  else:cb.set_ticks([-limits[j],-.1,0,.1,limits[j]]);cb.set_ticklabels([f'{-limits[j]:.2g}','−0.1','0','0.1',f'{limits[j]:.2g}']);cb.set_label('ζ/f₀ (+ cyclonic, − anticyclonic)')
  axs[row,j].set(title=('Cap jet' if row<2 else 'No jet')+' · '+('tracer' if tracer else 'vorticity')+f' · z/H={heights[j]:.5f}',xlabel='x/H',ylabel='y/H')
title=fig.suptitle('',y=.982)
def update(t):
 for block,label in enumerate(['jet','nojet']):
  with h5py.File(frames[label][t]) as h:q=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
  assert np.isfinite(q).all() and np.isfinite(tr).all()
  for j in range(3):images[block*6+j].set_data(q[j]);images[block*6+3+j].set_data(tr[j])
 title.set_text(f'Anelastic heated β channel · Ω=3 · t={t} (nondimensional)\n128 × 256 × 80 · WENO9 / no closure · fixed scales · movie sampling Δt=5')
# Check that all axes titles and labels fit before encoding the complete movie.
for t in [0,1000]:
 update(t);fig.canvas.draw();renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
 for ax in fig.axes:
  for text in [ax.title,ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()]:
   if text.get_visible() and text.get_text():
    b=text.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=width and b.y1<=height,(text.get_text(),b.bounds,width,height)
 assert width%2==height%2==0
fig.savefig(out/'anelastic_layout_check.png',dpi=90)
writer=FFMpegWriter(fps=8,codec='libx264',bitrate=1500,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
movie=out/'anelastic_0_1000.mp4'
with writer.saving(fig,str(movie),dpi=90):
 for t in range(0,1001,5):update(t);writer.grab_frame()
fig.savefig(out/'anelastic_1000.png',dpi=90);plt.close(fig)
n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==201
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
for name in ['anelastic_0_1000.mp4','anelastic_1000.png']:
 p=out/name;m['assets'][name]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
m['anelastic_layout']='Explicit subplot margins; axis text bounds checked at t0 and1000;1260x1350 pixels';(out/'media.json').write_text(json.dumps(m,indent=2));print('ANELASTIC_LAYOUT_AND_DECODE_VERIFIED',flush=True)

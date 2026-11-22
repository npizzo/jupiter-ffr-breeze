import sys,json,subprocess
from pathlib import Path
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
root=Path('/work/pi_nicholas_pizzo_uri_edu/jupiter-ffr-breeze');src=root/'data/omega3-anelastic-v4-run-2113277/parent';out=root/sys.argv[1];out.mkdir(parents=True,exist_ok=False)
frames=sorted(src.glob('frame_*.jld2'));assert len(frames)==11
limits=np.zeros(3)
for p in frames:
 with h5py.File(p) as f:limits=np.maximum(limits,np.max(np.abs(f['vorticity_slices'][...]/float(f['f'][()])),axis=(1,2)))
limits=np.maximum(limits,.05)
fig,axes=plt.subplots(2,3,figsize=(13.6,8),layout='constrained');ims=[]
with h5py.File(frames[0]) as f:z=f['video_heights'][...]
for k in range(3):
 ims.append(axes[0,k].imshow(np.zeros((256,256)),origin='lower',extent=[0,4,0,4],cmap='RdBu_r',norm=AsinhNorm(linear_width=.05,vmin=-limits[k],vmax=limits[k])))
 fig.colorbar(ims[-1],ax=axes[0,k],label='ζ / f₀ (fixed asinh)');axes[0,k].set_title(f'z/H = {z[k]:.5f}')
for k in range(3):
 ims.append(axes[1,k].imshow(np.zeros((256,256)),origin='lower',extent=[0,4,0,4],cmap='BrBG',norm=AsinhNorm(linear_width=.01,vmin=0,vmax=1)))
 cb=fig.colorbar(ims[-1],ax=axes[1,k],label='Tracer concentration');cb.set_ticks([0,.01,.05,.1,.2,.5,1]);axes[1,k].set_title('Passive tracer; fixed asinh width 0.01')
for ax in axes.flat:ax.set_xlabel('x/H');ax.set_ylabel('y/H')
title=fig.suptitle('')
writer=FFMpegWriter(fps=2,codec='libx264',extra_args=['-crf','20','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(out/'pilot.mp4'),dpi=100):
 for p in frames:
  with h5py.File(p) as f:
   t=float(f['t'][()]);v=f['vorticity_slices'][...]/float(f['f'][()]);c=f['tracer_slices'][...]
  assert np.all(np.isfinite(v)) and np.all(np.isfinite(c))
  for k in range(3):ims[k].set_data(v[k]);ims[k+3].set_data(c[k])
  title.set_text(f'Ω = 3 anelastic validation pilot | t = {t:.2f} | Not an LSV formation test yet')
  writer.grab_frame()
 fig.savefig(out/'poster.png',dpi=100)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(out/'pilot.mp4'))
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(out/'pilot.mp4'),'-f','null','-'],check=True)
assert count==11 and abs(duration-5.5)<.1
(out/'media.json').write_text(json.dumps(dict(status='PASS',frames=11,duration=duration,heights=z.tolist(),vorticity_limits=limits.tolist(),vorticity_scale='asinh width 0.05, fixed per-depth symmetric full extrema over all frames',tracer_scale='asinh width 0.01; fixed range 0–1; BrBG',source=str(src.relative_to(root))),indent=2)+'\n')
print('MEDIA_PASS',str(out),flush=True)

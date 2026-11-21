"""Unity-only fixed-scale 3-depth movie; retains actual passive concentration."""
from pathlib import Path
import sys,json,subprocess
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm
import imageio_ffmpeg
source,out=map(Path,sys.argv[1:3]);out.mkdir(parents=True,exist_ok=True)
files=sorted(source.glob('frame_*.jld2'));assert len(files)>1
vmax=0.;cmin=float('inf');cmax=-float('inf');times=[]
for p in files:
 with h5py.File(p) as h:
  t=float(h['t'][()]);times.append(t);vmax=max(vmax,float(abs(h['vorticity_slices'][()]).max()));c=h['tracer_slices'][()];cmin=min(cmin,float(c.min()));cmax=max(cmax,float(c.max()))
assert all(np.diff(times)>0)
channel=times[0]>1e6;f=0.0002880957739764384 if channel else .5
extent=(0,6000,-3000,3000) if channel else (0,2,0,2)
center=.0793082 if channel else .1;width=.01
# Include the actual absolute value even when every stored tracer value is equal.
lo=min(cmin,center-.01);hi=max(cmax,center+.01)
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ffmpeg
fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained');ims=[]
with h5py.File(files[0]) as h:heights=h['video_heights'][()] if 'video_heights' in h else h['video_heights_requested_m'][()]
for j in range(3):
 ims.append(axs[0,j].imshow(np.zeros((8,8)),origin='lower',extent=extent,cmap='RdBu_r',norm=AsinhNorm(linear_width=.02,vmin=-vmax/f,vmax=vmax/f),interpolation='nearest'))
 axs[0,j].set_title(f'z={heights[j]/1000:.2f} km' if channel else f'z={heights[j]:.4f}')
for j in range(3):ims.append(axs[1,j].imshow(np.zeros((8,8)),origin='lower',extent=extent,cmap='BrBG',norm=AsinhNorm(linear_width=width,vmin=lo-center,vmax=hi-center),interpolation='nearest'))
fig.colorbar(ims[0],ax=axs[0,:].tolist(),label='ζ/f; positive = cyclonic')
cb=fig.colorbar(ims[3],ax=axs[1,:].tolist(),label=f'Absolute passive concentration; asinh width {width}')
ticks=np.linspace(lo,hi,5);cb.set_ticks(ticks-center);cb.set_ticklabels([f'{x:.4g}' for x in ticks])
for ax in axs.flat:ax.set_xlabel('x [km]' if channel else 'x');ax.set_ylabel('y [km]' if channel else 'y')
writer=FFMpegWriter(fps=10,codec='libx264',bitrate=1100,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart'])
movie=out/'evolution.mp4'
with writer.saving(fig,str(movie),dpi=90):
 for p,t in zip(files,times):
  with h5py.File(p) as h:
   z=h['vorticity_slices'][()]/f;c=h['tracer_slices'][()]
  for j in range(3):ims[j].set_data(z[j]);ims[j+3].set_data(c[j]-center)
  fig.suptitle(f'{source.parent.name} · day {t/86400:.3f}' if channel else f'{source.parent.name} · t={t:.2f}')
  writer.grab_frame()
subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-f','null','-'],check=True)
assert movie.stat().st_size<24_000_000
(out/'movie.json').write_text(json.dumps({'frame_count':len(files),'start':times[0],'end':times[-1],'heights':heights.tolist(),'f':f,'vorticity_limit':vmax/f,'tracer_center':center,'tracer_width':width,'tracer_data_range':[cmin,cmax],'tracer_color_range':[lo,hi],'fixed_across_depth_and_time':True,'fully_decoded':True},indent=2))
print('MOVIE_DECODE_PASS',movie,flush=True)

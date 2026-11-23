import json,hashlib,subprocess
from pathlib import Path
import h5py,numpy as np,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
root=Path('data/beta-shear-v1-pilot-2134379/parent')
out=Path('data/beta-shear-v1-postmortem');out.mkdir(exist_ok=True)
frames=sorted(root.glob('frame_*.jld2'));lim=np.zeros(3)
for p in frames:
 with h5py.File(p) as h:
  z=h['vorticity_slices'][:]/float(h['f'][()]);assert np.isfinite(z).all();lim=np.maximum(lim,np.max(abs(z),axis=(1,2)));heights=h['video_heights'][:]
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
fig,axs=plt.subplots(2,3,figsize=(10,10));fig.subplots_adjust(top=.89,bottom=.07,left=.06,right=.94,wspace=.48,hspace=.25)
ims=[]
for row in range(2):
 for j in range(3):
  norm=AsinhNorm(linear_width=.01,vmin=0,vmax=1) if row else AsinhNorm(linear_width=.02,vmin=-max(.05,lim[j]),vmax=max(.05,lim[j]))
  im=axs[row,j].imshow(np.zeros((256,128)),origin='lower',extent=[0,4,-4,4],cmap='BrBG' if row else 'RdBu_r',norm=norm,interpolation='nearest');ims.append(im)
  fig.colorbar(im,ax=axs[row,j],fraction=.06,pad=.04,label='Tracer' if row else 'ζ/f₀ (cyclonic +)')
  axs[row,j].set_title(f'z/H={heights[j]:.3f}');axs[row,j].set_xlabel('x/H');axs[row,j].set_ylabel('y/H')
title=fig.suptitle('');movie=out/'failed_jet_pilot_0_158.mp4'
w=FFMpegWriter(fps=12,codec='libx264',bitrate=650,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
with w.saving(fig,str(movie),dpi=80):
 for p in frames:
  with h5py.File(p) as h:
   t=float(h['t'][()]);z=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
  for j in range(3):ims[j].set_data(z[j]);ims[3+j].set_data(tr[j])
  title.set_text(f'FAILED unheated shear validation · t={t:g} (nondimensional)\n'+('Before Mach acceptance limit crossed' if t<51.5 else 'REJECTED STATE: Mach acceptance limit crossed at t=51.5')+'\nFull vorticity includes imposed jets; fixed scales throughout')
  w.grab_frame()
fig.savefig(out/'failed_jet_pilot_final.png',dpi=80);plt.close(fig)
n,dur=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==len(frames)
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
(out/'media.json').write_text(json.dumps(dict(status='PASS',simulation_status='FAILED',frames=n,duration_seconds=dur,bytes=movie.stat().st_size,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),vorticity_limits=lim.tolist(),vorticity_asinh_width=.02,tracer_asinh_width=.01,tracer_range=[0,1],heights=heights.tolist(),source=str(root)),indent=2))
print('MEDIA_PASS',n,dur,movie.stat().st_size)

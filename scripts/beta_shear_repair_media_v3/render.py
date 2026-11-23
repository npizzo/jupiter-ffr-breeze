"""Render existing three-depth fields on Unity only; no model extrapolation."""
import sys,json,hashlib,subprocess,tomllib
from pathlib import Path
import h5py,numpy as np,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.ticker import ScalarFormatter
from matplotlib.animation import FFMpegWriter
out=Path(sys.argv[1]);out.mkdir(exist_ok=True,parents=True)
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe();media=[]
for source in map(Path,sys.argv[2:]):
 for root in sorted(source.iterdir()):
  if root.name not in ['weno','centered_momentum','centered_both','cap_anticyclonic','cap_cyclonic','traditional_anticyclonic','traditional_cyclonic']:continue
  meta=tomllib.loads((root/'metadata.toml').read_text())
  frames=sorted(root.glob('frame_*.jld2'))
  if len(frames)<2:continue
  lim=np.zeros(3)
  for p in frames:
   with h5py.File(p) as h:
    z=h['vorticity_slices'][:]/float(h['f'][()]);assert np.isfinite(z).all();lim=np.maximum(lim,np.max(abs(z),axis=(1,2)));heights=h['video_heights'][:];lx=float(h['lx'][()]);ly=float(h['ly'][()])
  fig,axs=plt.subplots(2,3,figsize=(12,7.5),layout='constrained');ims=[]
  for row in range(2):
   for j in range(3):
    norm=AsinhNorm(linear_width=.01,vmin=0,vmax=1) if row else AsinhNorm(linear_width=.002,vmin=-max(.0001,lim[j]),vmax=max(.0001,lim[j]))
    im=axs[row,j].imshow(np.zeros((2,2)),origin='lower',extent=[0,lx,-ly/2,ly/2],aspect='auto',cmap='BrBG' if row else 'RdBu_r',norm=norm,interpolation='nearest');ims.append(im)
    cb=fig.colorbar(im,ax=axs[row,j],fraction=.05,pad=.02,label='Tracer' if row else 'ζ/f₀ (cyclonic +)')
    if row:cb.set_ticks([0,.01,.1,1],labels=['0','0.01','0.1','1'])
    else:
     bound=max(.0001,lim[j]);cb.set_ticks([-bound,-bound/2,0,bound/2,bound]);formatter=ScalarFormatter(useMathText=True);formatter.set_powerlimits((-2,2));cb.formatter=formatter;cb.update_ticks()
    axs[row,j].set(title=f'z/H={heights[j]:.3f}',xlabel='x/H',ylabel='y/H')
  title=fig.suptitle('');movie=out/(root.name+'.mp4');w=FFMpegWriter(fps=4,codec='libx264',bitrate=600,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
  with w.saving(fig,str(movie),dpi=80):
   for p in frames:
    with h5py.File(p) as h:t=float(h['t'][()]);z=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
    assert np.isfinite(tr).all()
    for j in range(3):ims[j].set_data(z[j]);ims[3+j].set_data(tr[j])
    title.set_text(f'{root.name.replace("_"," ")} · t={t:g} · {meta["status"]} at t={meta["final_time"]:.3f}\nBounded validation; no bottom heating. Fixed scales; axes stretched for visibility.')
    w.grab_frame()
  fig.savefig(out/(root.name+'.png'),dpi=80);plt.close(fig)
  n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==len(frames)
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
  media.append(dict(simulation_status=meta["status"],simulation_final_time=meta["final_time"],movie=movie.name,source=str(root),frames=n,duration_seconds=duration,final_time=t,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),bytes=movie.stat().st_size,heights=heights.tolist(),vorticity_limits=np.maximum(lim,.0001).tolist(),vorticity_asinh_width=.002,tracer_range=[0,1],tracer_asinh_width=.01))
  print('MEDIA_PASS',root.name,n,t,flush=True)
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n')

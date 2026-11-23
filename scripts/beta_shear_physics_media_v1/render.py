from pathlib import Path
import json,hashlib,subprocess,shutil
import h5py,numpy as np,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter
root=Path('data/beta-shear-physics-v1-run-2162603');out=Path('rendered/beta-shear-physics100');out.mkdir(parents=True,exist_ok=False)
frames={};limits=np.zeros(3);trmin=1.;trmax=0.
for case in ['jet','nojet']:
 frames[case]={}
 for suffix in ['pilot','to100']:
  for p in sorted((root/(case+'_'+suffix)).glob('frame_*.jld2')):
   with h5py.File(p) as h:
    t=float(h['t'][()]);assert abs(t-round(t))<1e-7
    if round(t) in frames[case]:
     with h5py.File(frames[case][round(t)]) as old:
      assert np.array_equal(old['vorticity_slices'][:],h['vorticity_slices'][:]);assert np.array_equal(old['tracer_slices'][:],h['tracer_slices'][:])
    frames[case][round(t)]=p
    q=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
    assert np.isfinite(q).all() and np.isfinite(tr).all()
    limits=np.maximum(limits,np.max(abs(q),axis=(1,2)));trmin=min(trmin,float(tr.min()));trmax=max(trmax,float(tr.max()))
    heights=h['video_heights'][:];lx=float(h['lx'][()]);ly=float(h['ly'][()])
 assert sorted(frames[case])==list(range(101))
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
fig,axs=plt.subplots(4,3,figsize=(14,15),layout='constrained');images=[]
center=.1;width=.005
normtr=FuncNorm((lambda x:np.arcsinh((x-center)/width),lambda x:center+width*np.sinh(x)),vmin=0,vmax=1)
for row in range(4):
 for j in range(3):
  tracer=row%2==1;norm=normtr if tracer else AsinhNorm(linear_width=.01,vmin=-limits[j],vmax=limits[j])
  im=axs[row,j].imshow(np.zeros((64,32)),origin='lower',extent=[0,lx,-ly/2,ly/2],aspect='auto',interpolation='nearest',cmap='BrBG' if tracer else 'RdBu_r',norm=norm);images.append(im)
  cb=fig.colorbar(im,ax=axs[row,j],fraction=.045,pad=.02)
  if tracer:cb.set_ticks([0,.05,.09,.1,.11,.15,1],labels=['0','.05','.09','.10','.11','.15','1']);cb.set_label('Tracer concentration C')
  else:
   ticks=[-limits[j],-.1,0,.1,limits[j]] if limits[j]>.11 else np.linspace(-limits[j],limits[j],5)
   cb.set_ticks(ticks);cb.set_ticklabels([f'{v:.3g}' for v in ticks]);cb.set_label('ζ/f₀: cyclonic +, anticyclonic −')
  label='Modest cap jet' if row<2 else 'No jet'
  axs[row,j].set(title=f'{label} · '+('tracer' if tracer else 'vorticity')+f' · z/H={heights[j]:.5f}',xlabel='x/H',ylabel='y/H')
title=fig.suptitle('');writer=FFMpegWriter(fps=8,codec='libx264',bitrate=1000,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
movie=out/'jet_nojet_0_100.mp4'
with writer.saving(fig,str(movie),dpi=90):
 for t in range(101):
  for block,case in enumerate(['jet','nojet']):
   with h5py.File(frames[case][t]) as h:q=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
   for j in range(3):images[block*6+j].set_data(q[j]);images[block*6+3+j].set_data(tr[j])
  title.set_text(f'Heated compressible beta channel · Ω=3 · t={t:g} (nondimensional)\nMatched eddies; fixed scales across cases and time · coarse exploratory grid 32 × 64 × 80')
  writer.grab_frame()
  if t in [0,50,100]:fig.savefig(out/f'frame_{t:03d}.png',dpi=90)
plt.close(fig)
n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==101
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
summary=json.loads((root/'PHYSICS_REVIEW.json').read_text())
# First resumed output is one step after t20, not a same-time duplicate; relabel that difference accurately.
for c in summary['cases'].values():c['first_resumed_sample_increment']=c.pop('restart_budget_difference')
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
meta={'status':'MEDIA_VERIFIED','source_job':2162603,'source':str(root),'simulation_times':list(range(101)),'frames':n,'fps':8,'duration_seconds':duration,'full_decode_passed':True,'heights':heights.tolist(),'layout':'4x3: jet vorticity,tracer; nojet vorticity,tracer','axes':'x/H0..4,y/H-4..4; panels stretched for readability','vorticity_normalization':'zeta/f0; total vorticity including jet; positive cyclonic in Northern Hemisphere','vorticity_limits_per_depth':limits.tolist(),'vorticity_asinh_width':.01,'tracer_transform':'asinh((C-0.1)/0.005), fixed across all cases, depths and times; ticks show absolute concentration','tracer_range':[0,1],'tracer_saved_extrema':[trmin,trmax],'timeline':'All101 integer states; pilot0..20 plus continuation21..100; no temporal interpolation','assets':{}}
for p in out.iterdir():
 if p.is_file():meta['assets'][p.name]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(out/'media.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))

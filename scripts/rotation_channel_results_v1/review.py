"""Read-only channel-pair review and rendering; execute on Unity."""
import sys,json,csv,tomllib,hashlib,subprocess
from pathlib import Path
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
roots={'beta':Path('data/rotation-channel-v1-beta-2109125/parent'),'control':Path('data/rotation-channel-v1-control-2109124/parent')}
frames={};summaries={};limits=np.zeros(3);records={}
for case,root in roots.items():
 meta=tomllib.loads((root/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==1000
 ff={};rr=[]
 for p in sorted(root.glob('frame_*.jld2')):
  with h5py.File(p) as h:
   t=float(h['t'][()]);z=h['vorticity_slices'][()].astype(float);tr=h['tracer_slices'][()];f=float(h['f'][()]);heights=h['video_heights'][()].tolist()
  assert z.shape==(3,256,256) and np.isfinite(z).all() and np.isfinite(tr).all()
  assert abs(t-round(t))<1e-8 and t not in ff
  ff[t]=p;limits=np.maximum(limits,np.max(abs(z/f),axis=(1,2)))
  r={'time':t}
  for j in range(3):
   q=z[j];r[f'skew_{j}']=float(np.mean(q**3)/max(np.mean(q*q)**1.5,1e-99));r[f'anticyclonic_enstrophy_share_{j}']=float(np.sum(q[q<0]**2)/max(np.sum(q*q),1e-99))
  rr.append(r)
 assert np.array_equal(sorted(ff),np.arange(1001)),(case,len(ff),min(ff),max(ff))
 frames[case]=ff;records[case]=rr
 with (root/'diagnostics.csv').open() as h:dr=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(h)]
 assert all(all(np.isfinite(v) for v in r.values()) for r in dr)
 with (out/f'{case}_planes.csv').open('w') as h:w=csv.DictWriter(h,fieldnames=rr[0]);w.writeheader();w.writerows(rr)
 summaries[case]={'source':str(root),'final_time':1000,'frames':len(ff),'f0':f,'beta':meta.get('beta',0),'actual_heights':heights,'late_window_950_1000':{k:float(np.mean([r[k] for r in rr if r['time']>=950])) for k in rr[0] if k!='time'},'final_diagnostics':dr[-1],'max_relative_mass_change':max(abs(r['mass']/dr[0]['mass']-1) for r in dr)}
 print(case,'DATA_PASS',json.dumps(summaries[case]),flush=True)
summary={'cases':summaries,'geometry':'x periodic; y impermeable, stress-free, insulated walls; 4 x 4 x 1','omega':3,'latitude_degrees_north':67.5,'depth_mapping_km':100,'width_mapping_km':400,'beta_fractional_f_change_across_box':summaries['beta']['beta']*4/f,'limitations':['One realization per case.','Plane enstrophy and skewness are not coherent-vortex tracking or inverse-cascade diagnostics.','Channel boundaries and wall-compatible initial eddies differ from the original doubly periodic runs.','Original compressible energy-budget residual remains unresolved.','Nondimensional times, not seconds or Earth days.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
fig,axs=plt.subplots(4,3,figsize=(12,12));fig.subplots_adjust(left=.055,right=.94,bottom=.05,top=.93,wspace=.32,hspace=.45);ims=[];limits=np.maximum(limits,.05)
for row in range(4):
 case='beta' if row<2 else 'control';tracer=row%2==1
 for j in range(3):
  norm=AsinhNorm(linear_width=.01,vmin=0,vmax=1) if tracer else AsinhNorm(linear_width=.02,vmin=-limits[j],vmax=limits[j])
  im=axs[row,j].imshow(np.zeros((256,256)),origin='lower',extent=[0,4,-2,2],interpolation='nearest',cmap='BrBG' if tracer else 'RdBu_r',norm=norm);ims.append(im)
  label='Tracer' if tracer else 'ζ/f₀ (cyclonic +)';cb=fig.colorbar(im,ax=axs[row,j],fraction=.046,pad=.025,label=label)
  if tracer:cb.set_ticks([0,.01,.05,.1,.2,.5,1]);cb.set_ticklabels(['0','.01','.05','.1','.2','.5','1'])
  axs[row,j].set_title(f'{"Beta" if case=="beta" else "f-plane control"} · z/H={heights[j]:.3f}',fontsize=10)
  axs[row,j].set_xlabel('x/H');axs[row,j].set_ylabel('y/H')
title=fig.suptitle('');movie=out/'channel_pair_0_1000.mp4'
writer=FFMpegWriter(fps=20,codec='libx264',bitrate=1100,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=90):
 for t in range(1001):
  for c,case in enumerate(['beta','control']):
   with h5py.File(frames[case][t]) as h:z=h['vorticity_slices'][()]/f;tr=h['tracer_slices'][()]
   for j in range(3):ims[c*6+j].set_data(z[j]);ims[c*6+3+j].set_data(tr[j])
  title.set_text(f'Ω=3 channel comparison · t={t} (nondimensional)\nBeta above; matching f-plane control below · shared fixed color scales')
  writer.grab_frame()
  if t%100==0:print('RENDER',t,flush=True)
fig.savefig(out/'channel_pair_1000.png',dpi=100);plt.close(fig)
n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==1001 and abs(duration-50.05)<.1
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
media={'status':'PASS','frames':n,'duration_seconds':duration,'sha256':hashlib.sha256(movie.read_bytes()).hexdigest(),'bytes':movie.stat().st_size,'vorticity_limits_per_depth':limits.tolist(),'vorticity_asinh_width':.02,'tracer_asinh_width':.01,'tracer_range':[0,1],'shared_scales_between_cases':True,'location':'Unity'}
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n');print('PUBLICATION_READY',json.dumps(media),flush=True)

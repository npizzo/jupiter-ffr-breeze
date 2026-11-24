from pathlib import Path
import csv,json,hashlib,subprocess,tomllib
import h5py,numpy as np,imageio_ffmpeg
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter
root=Path.cwd();out=root/'rendered/weno9-lsv1000-nojet-v2';out.mkdir(parents=True,exist_ok=False)
dirs=[root/'data/weno9-startup-v2-pilot-2177387_0/parent',root/'data/weno9-lsv100-v3-run-2181508_0/parent',root/'data/weno9-lsv1000-v4-run-2187099_0/parent']
meta=tomllib.loads((dirs[-1]/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==1000
frames={};ds={};ss={};limits=np.zeros(3)
for d in dirs:
 for filename,target in [('diagnostics.csv',ds),('mach_spatial.csv',ss)]:
  for r in csv.DictReader((d/filename).open()):target[float(r['time'])]={k:float(v) for k,v in r.items()}
 for p in sorted(d.glob('frame_*.jld2')):
  with h5py.File(p) as h:
   t=float(h['t'][()]);assert abs(t-round(t))<1e-7
   if round(t) in frames:
    with h5py.File(frames[round(t)]) as old:assert np.array_equal(old['vorticity_slices'][:],h['vorticity_slices'][:])
   frames[round(t)]=p;q=h['vorticity_slices'][:]/float(h['f'][()]);assert np.isfinite(q).all();limits=np.maximum(limits,np.max(abs(q),axis=(1,2)));heights=h['video_heights'][:];lx=float(h['lx'][()]);ly=float(h['ly'][()])
assert sorted(frames)==list(range(1001));series=[];smooth_final=None
for t,p in sorted(frames.items()):
 with h5py.File(p) as h:
  u=h['ubar'][:];v=h['vbar'][:];f=float(h['f'][()]);ny,nx=u.shape;dx=lx/nx;dy=ly/ny
  un=u-u.mean(axis=1,keepdims=True);vn=v-v.mean(axis=1,keepdims=True)
  zeta=((np.roll(vn,-1,axis=1)-np.roll(vn,1,axis=1))/(2*dx)-np.gradient(un,dy,axis=0))/f
  smooth=gaussian_filter(zeta,(.25/dy,.25/dx),mode=('reflect','wrap'));y=(np.arange(ny)+.5)*dy-ly/2;inside=abs(y)<3.5;a=smooth[inside];en=np.sum(a*a);q=h['vorticity_slices'][:]/f
  row={'time':t,'depthmean_total_ke_per_mass':float(.5*np.mean(u*u+v*v)),'depthmean_nonzonal_ke_per_mass':float(.5*np.mean(un*un+vn*vn)),'smoothed_interior_negative_enstrophy_share':float(np.sum(a*a*(a<0))/en),'smoothed_interior_min_zeta_over_f':float(a.min()),'smoothed_interior_max_zeta_over_f':float(a.max())}
  for k in range(3):
   b=q[k]-q[k].mean(axis=1,keepdims=True);b=b[inside];e=np.sum(b*b);row[f'depth{k}_nonzonal_negative_enstrophy_share']=float(np.sum(b*b*(b<0))/e)
  series.append(row)
  if t==1000:smooth_final=smooth
late=[r for r in series if r['time']>=900]
summary={'status':'COMPRESSIBLE_NOJET_COMPLETE_1000','time_units':'nondimensional','source_jobs':[2177387,2181508,2187099],'final_diagnostics':ds[max(ds)],'final_spatial_mach':ss[max(ss)],'late_900_1000':{k:float(np.mean([r[k] for r in late])) for k in late[0] if k!='time'},'milestones':[r for r in series if r['time'] in [0,100,250,500,750,1000]],'configuration':{k:meta[k] for k in ['nx','ny','nz','equations','jet_speed','initial_rms','advection_order','closure']},'diagnostic_notes':['No imposed jet and no jet momentum relaxation in this case.','Spatial diagnostic removes zonal mean then Gaussian-smooths depth-mean vorticity with sigma0.25H, excluding |y|>=3.5H; negative share is not a vortex count or flux.','Initial conditions already contain large mixed-sign eddies.','Temperature-profile drift and compressible energy residual are reported; completion is not quantitative validation.','All1001 integer states analyzed; movie shows every5 time units,201 frames; no temporal interpolation.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (out/'large_scale_metrics.csv').open('w') as io:
 w=csv.DictWriter(io,fieldnames=list(series[0]));w.writeheader();w.writerows(series)
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
fig,axs=plt.subplots(2,3,figsize=(13.2,8.4),layout='constrained');images=[]
normtr=FuncNorm((lambda x:np.arcsinh((x-.1)/.005),lambda x:.1+.005*np.sinh(x)),vmin=0,vmax=1)
for row in range(2):
 for j in range(3):
  tracer=row==1;norm=normtr if tracer else AsinhNorm(linear_width=.01,vmin=-limits[j],vmax=limits[j])
  im=axs[row,j].imshow(np.zeros((256,128)),origin='lower',extent=[0,lx,-ly/2,ly/2],aspect='auto',interpolation='nearest',cmap='BrBG' if tracer else 'RdBu_r',norm=norm);images.append(im);cb=fig.colorbar(im,ax=axs[row,j],fraction=.045,pad=.02)
  if tracer:cb.set_ticks([0,.05,.09,.1,.11,.15,1],labels=['0','.05','.09','.10','.11','.15','1']);cb.set_label('Tracer C')
  else:cb.set_ticks([-limits[j],-.1,0,.1,limits[j]]);cb.set_ticklabels([f'{-limits[j]:.2g}','−0.1','0','0.1',f'{limits[j]:.2g}']);cb.set_label('ζ/f₀ (+ cyclonic, − anticyclonic)')
  axs[row,j].set(title=('Tracer' if tracer else 'Vorticity')+f' · z/H={heights[j]:.5f}',xlabel='x/H',ylabel='y/H')
assert all(round(x*90)%2==0 for x in fig.get_size_inches())
title=fig.suptitle('');writer=FFMpegWriter(fps=8,codec='libx264',bitrate=1400,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart']);movie=out/'compressible_nojet_0_1000.mp4'
with writer.saving(fig,str(movie),dpi=90):
 for t in range(0,1001,5):
  with h5py.File(frames[t]) as h:q=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
  assert np.isfinite(tr).all()
  for j in range(3):images[j].set_data(q[j]);images[3+j].set_data(tr[j])
  title.set_text(f'Compressible heated β channel · no imposed jet · Ω=3 · t={t} (nondimensional)\n128 × 256 × 80 · WENO9 / no closure · fixed scales · movie sampling Δt=5')
  writer.grab_frame()
  if t in [0,100,500,1000]:fig.savefig(out/f'frame_{t:04d}.png',dpi=90)
plt.close(fig)
n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==201
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
fig,ax=plt.subplots(figsize=(5,8),layout='constrained');lim=float(abs(smooth_final).max());im=ax.imshow(smooth_final,origin='lower',extent=[0,4,-4,4],cmap='RdBu_r',norm=AsinhNorm(linear_width=.01,vmin=-lim,vmax=lim));cb=fig.colorbar(im,ax=ax,label='Smoothed nonzonal depth-mean ζ/f₀');cb.set_ticks([-lim,-.01,0,.01,lim]);cb.set_ticklabels([f'{-lim:.3f}','−0.01','0','0.01',f'{lim:.3f}']);ax.set(title='Compressible, no jet · t = 1000\nZonal mean removed; Gaussian σ = 0.25H',xlabel='x/H',ylabel='y/H');fig.savefig(out/'depthmean_1000.png',dpi=110);plt.close(fig)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
m={'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames':201,'fps':8,'duration':duration,'times':list(range(0,1001,5)),'heights':heights.tolist(),'vorticity_limits':limits.tolist(),'vorticity_asinh_width':.01,'tracer_transform':'asinh((C-0.1)/0.005),range0..1;fixed across time/depth','assets':assets};(out/'media.json').write_text(json.dumps(m,indent=2));print(json.dumps(summary),flush=True)

from pathlib import Path
import csv,json,hashlib,subprocess,tomllib
import h5py,numpy as np,imageio_ffmpeg
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter
root=Path.cwd();out=root/'rendered/weno9-lsv1000';out.mkdir(parents=True,exist_ok=False)
frames={};limits=np.zeros(3);results={};series=[]
labels=['compressible_nojet','compressible_jet','anelastic_nojet','anelastic_jet']
for i,label in enumerate(labels):
 dirs=[root/f'data/weno9-startup-v2-pilot-2177387_{i}/parent',root/f'data/weno9-lsv100-v3-run-2181508_{i}/parent',root/f'data/weno9-lsv1000-v4-run-2187099_{i}/parent']
 meta=tomllib.loads((dirs[-1]/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==1000
 frames[label]={};ds={};ss={}
 for d in dirs:
  for filename,target in [('diagnostics.csv',ds),('mach_spatial.csv',ss)]:
   for r in csv.DictReader((d/filename).open()):target[float(r['time'])]={k:float(v) for k,v in r.items()}
  for p in sorted(d.glob('frame_*.jld2')):
   with h5py.File(p) as h:
    t=float(h['t'][()]);assert abs(t-round(t))<1e-7
    if int(round(t)) in frames[label]:
     with h5py.File(frames[label][int(round(t))]) as old:assert np.array_equal(old['vorticity_slices'][:],h['vorticity_slices'][:])
    frames[label][int(round(t))]=p
    q=h['vorticity_slices'][:]/float(h['f'][()]);assert np.isfinite(q).all();limits=np.maximum(limits,np.max(abs(q),axis=(1,2)))
    heights=h['video_heights'][:];lx=float(h['lx'][()]);ly=float(h['ly'][()])
 assert sorted(frames[label])==list(range(1001))
 for t,p in frames[label].items():
  with h5py.File(p) as h:
   u=h['ubar'][:];v=h['vbar'][:];f=float(h['f'][()]);ny,nx=u.shape;dx=lx/nx;dy=ly/ny
   # A spatial diagnostic, not a spectral flux or inverse-cascade measurement.
   un=u-u.mean(axis=1,keepdims=True);vn=v-v.mean(axis=1,keepdims=True)
   zeta=((np.roll(vn,-1,axis=1)-np.roll(vn,1,axis=1))/(2*dx)-np.gradient(un,dy,axis=0))/f
   smooth=gaussian_filter(zeta,(.25/dy,.25/dx),mode=('reflect','wrap'))
   y=(np.arange(ny)+.5)*dy-ly/2;interior=abs(y)<3.5;a=smooth[interior];en=np.sum(a*a)
   q=h['vorticity_slices'][:]/f
   row={'case':label,'time':t,'depthmean_total_ke_per_mass':float(.5*np.mean(u*u+v*v)),'depthmean_nonzonal_ke_per_mass':float(.5*np.mean(un*un+vn*vn)),'smoothed_interior_negative_enstrophy_share':float(np.sum(a*a*(a<0))/en) if en else 0.,'smoothed_interior_min_zeta_over_f':float(a.min()),'smoothed_interior_max_zeta_over_f':float(a.max())}
   for k in range(3):
    b=q[k]-q[k].mean(axis=1,keepdims=True);b=b[interior];e=np.sum(b*b)
    row[f'depth{k}_nonzonal_negative_enstrophy_share']=float(np.sum(b*b*(b<0))/e)
   series.append(row)
 late=[r for r in series if r['case']==label and r['time']>=900]
 results[label]={'final_diagnostics':ds[max(ds)],'final_spatial_mach':ss[max(ss)],'late_900_1000':{k:float(np.mean([r[k] for r in late])) for k in late[0] if k not in ['case','time']},'initial_large_scale':[r for r in series if r['case']==label and r['time']==0][0],'final_large_scale':[r for r in series if r['case']==label and r['time']==1000][0],'metadata':{k:meta[k] for k in ['nx','ny','nz','equations','jet_speed','initial_rms','advection_order','closure']}}
summary={'status':'ALL_FOUR_COMPLETE_1000','time_units':'nondimensional','source_jobs':[2177387,2181508,2187099],'cases':results,'diagnostic_notes':['Depth-mean uses the same reference-density weighting in both equations.','Large-scale diagnostic removes the zonal mean then Gaussian-smooths vorticity with sigma=0.25H; excludes |y|>=3.5H. This is an explicitly chosen spatial scale, not a spectral flux measurement.','Negative zeta/f is anticyclonic. Negative-enstrophy share alone does not establish a coherent LSV.','Compressible energy residual is recorded; anelastic energy_proxy is NOT a conserved anelastic energy invariant, so residuals cannot be compared as losses.','No proof of inverse transfer is inferred from morphology or relative energy fractions.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (out/'large_scale_metrics.csv').open('w') as io:
 w=csv.DictWriter(io,fieldnames=list(series[0]));w.writeheader();w.writerows(series)
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':9})
center=.1;width=.005
normtr=FuncNorm((lambda x:np.arcsinh((x-center)/width),lambda x:center+width*np.sinh(x)),vmin=0,vmax=1)
for eq in ['compressible','anelastic']:
 fig,axs=plt.subplots(4,3,figsize=(13.2,14),layout='constrained');images=[]
 for row in range(4):
  for j in range(3):
   tracer=row%2==1;norm=normtr if tracer else AsinhNorm(linear_width=.01,vmin=-limits[j],vmax=limits[j])
   im=axs[row,j].imshow(np.zeros((256,128)),origin='lower',extent=[0,lx,-ly/2,ly/2],aspect='auto',interpolation='nearest',cmap='BrBG' if tracer else 'RdBu_r',norm=norm);images.append(im)
   cb=fig.colorbar(im,ax=axs[row,j],fraction=.045,pad=.02)
   if tracer:cb.set_ticks([0,.05,.09,.1,.11,.15,1],labels=['0','.05','.09','.10','.11','.15','1']);cb.set_label('Tracer C')
   else:cb.set_label('ζ/f₀ (+ cyclonic, − anticyclonic)')
   axs[row,j].set(title=('Cap jet' if row<2 else 'No jet')+' · '+('tracer' if tracer else 'vorticity')+f' · z/H={heights[j]:.5f}',xlabel='x/H',ylabel='y/H')
 title=fig.suptitle(f'{eq.capitalize()} heated β channel · Ω=3 · t=1000 (nondimensional)\n128 × 256 × 80 · WENO9 / no closure · fixed scales · movie sampling Δt=5')
 fig.set_dpi(90);fig.canvas.draw();fig.set_layout_engine(None)
 assert all(round(x*90)%2==0 for x in fig.get_size_inches())
 writer=FFMpegWriter(fps=8,codec='libx264',bitrate=1500,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart']);movie=out/f'{eq}_0_1000.mp4'
 with writer.saving(fig,str(movie),dpi=90):
  for t in range(0,1001,5):
   for block,label in enumerate([eq+'_jet',eq+'_nojet']):
    with h5py.File(frames[label][t]) as h:q=h['vorticity_slices'][:]/float(h['f'][()]);tr=h['tracer_slices'][:]
    assert np.isfinite(tr).all()
    for j in range(3):images[block*6+j].set_data(q[j]);images[block*6+3+j].set_data(tr[j])
   title.set_text(f'{eq.capitalize()} heated β channel · Ω=3 · t={t} (nondimensional)\n128 × 256 × 80 · WENO9 / no closure · fixed scales · movie sampling Δt=5')
   writer.grab_frame()
   if t==1000:fig.savefig(out/f'{eq}_1000.png',dpi=90)
 plt.close(fig)
 n,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==201
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
 print(eq,'201 frames decoded',flush=True)
fields=[]
for label in labels:
 with h5py.File(frames[label][1000]) as h:
  u=h['ubar'][:];v=h['vbar'][:];f=float(h['f'][()]);ny,nx=u.shape;dx=lx/nx;dy=ly/ny
  u=u-u.mean(axis=1,keepdims=True);v=v-v.mean(axis=1,keepdims=True)
  q=((np.roll(v,-1,axis=1)-np.roll(v,1,axis=1))/(2*dx)-np.gradient(u,dy,axis=0))/f
  fields.append(gaussian_filter(q,(.25/dy,.25/dx),mode=('reflect','wrap')))
lim=max(float(abs(q).max()) for q in fields)
fig,axs=plt.subplots(2,2,figsize=(9,11),layout='constrained')
for ax,label,q in zip(axs.flat,labels,fields):
 im=ax.imshow(q,origin='lower',extent=[0,lx,-ly/2,ly/2],cmap='RdBu_r',norm=AsinhNorm(linear_width=.01,vmin=-lim,vmax=lim))
 ax.set(title=label.replace('_',' '),xlabel='x/H',ylabel='y/H');cb=fig.colorbar(im,ax=ax,label='Smoothed nonzonal depth-mean ζ/f₀');cb.set_ticks([-lim,-.01,0,.01,lim]);cb.set_ticklabels([f'{-lim:.3f}','−0.01','0','0.01',f'{lim:.3f}'])
fig.suptitle('t = 1000 · reference-density-weighted depth mean\nZonal mean removed; Gaussian σ = 0.25H; full extrema, no clipping')
fig.savefig(out/'depthmean_1000.png',dpi=110);plt.close(fig)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
meta={'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames_per_movie':201,'movie_sampling_interval':5,'integer_states_analyzed_per_case':1001,'fps':8,'times':[0,1000],'heights':heights.tolist(),'layout':'4x3 per movie: jet vorticity,tracer; nojet vorticity,tracer','vorticity_limits':limits.tolist(),'vorticity_asinh_width':.01,'tracer_transform':'asinh((C-0.1)/0.005), range0..1; fixed across cases,depths,time','assets':assets}
(out/'media.json').write_text(json.dumps(meta,indent=2));print(json.dumps(summary),flush=True)

from pathlib import Path
import json,csv,hashlib,subprocess
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter
root=Path.cwd();out=root/'rendered/weno9-fourstack';out.mkdir(parents=True,exist_ok=False)
old=json.loads((root/'rendered/weno9-lsv2000/media.json').read_text());assert old['full_decode_passed']
times=list(range(0,2001,5));labels=['Compressible · no jet','Compressible · cap jet','Anelastic · no jet','Anelastic · cap jet']
keys=['compressible_nojet','compressible_jet','anelastic_nojet','anelastic_jet'];Q=[];C=[];rows=[];profiles=[];tmin=np.inf;tmax=-np.inf
for i,key in enumerate(keys):
 paths={};diagnostics={}
 for job,prefix in [(2177387,'weno9-startup-v2-pilot'),(2181508,'weno9-lsv100-v3-run'),(2187099,'weno9-lsv1000-v4-run'),(2197612,'weno9-lsv2000-v5-run')]:
  d=root/f'data/{prefix}-{job}_{i}/parent'
  for r in csv.DictReader((d/'diagnostics.csv').open()):diagnostics[float(r['time'])]={k:float(v) for k,v in r.items()}
  for f in d.glob('frame_*.jld2'):
   t=round(int(f.stem[6:])/1000)
   if t in times:paths[t]=f
 assert sorted(paths)==times
 tmin=min(tmin,min(r['tracer_min'] for r in diagnostics.values()));tmax=max(tmax,max(r['tracer_max'] for r in diagnostics.values()))
 qi=np.empty((len(times),3,256,128),np.float32);ci=np.empty_like(qi);rho_initial=None
 for n,t in enumerate(times):
  with h5py.File(paths[t]) as h:
   qi[n]=h['vorticity_slices'][:]/float(h['f'][()]);ci[n]=h['tracer_slices'][:];z=h['z'][:];heights=h['video_heights'][:]
   u=h['ubar'][:];v=h['vbar'][:];mean_u=h['mean_u'][:];rho=h['mean_density'][:];T=h['mean_temperature'][:]
   if rho_initial is None:rho_initial=rho.copy()
   um=u.mean(axis=1,keepdims=True);vm=v.mean(axis=1,keepdims=True);ke=.5*np.mean(u*u+v*v);uke=.5*np.mean(um*um)
   r=dict(case=key,time=t,depthmean_total_ke=float(ke),depthmean_zonal_u_ke=float(uke),depthmean_zonal_u_fraction=float(uke/ke),depthmean_nonzonal_ke=float(.5*np.mean((u-um)**2+(v-vm)**2)),mean_density_profile_max_relative_change=float(np.max(abs(rho/rho_initial-1))))
   for j,zz in enumerate(heights):
    k=int(np.argmin(abs(z-zz)));a=ci[n,j];q=qi[n,j]
    r.update({f'depth{j}_mean_u_rms':float(np.sqrt(np.mean(mean_u[k]**2))),f'depth{j}_tracer_mean':float(a.mean(dtype=np.float64)),f'depth{j}_tracer_std':float(a.std(dtype=np.float64)),f'depth{j}_zonal_enstrophy_fraction':float(np.mean(q.mean(axis=1)**2)/np.mean(q*q))})
   rows.append(r)
   if t in [0,100,500,1000,1500,2000]:profiles.append(dict(case=key,time=t,z=z.tolist(),mean_density=rho.tolist(),mean_temperature=T.tolist(),heights=heights.tolist(),mean_u_at_sampled_heights=mean_u[[int(np.argmin(abs(z-zz))) for zz in heights]].tolist()))
 assert np.isfinite(qi).all() and np.isfinite(ci).all();Q.append(qi);C.append(ci)
 print(key,'401 FRAMES CACHED',flush=True)
with (out/'mean_flow_tracer.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=dict(cases={key:{'initial':next(r for r in rows if r['case']==key and r['time']==0),'final':next(r for r in rows if r['case']==key and r['time']==2000),'at1000':next(r for r in rows if r['case']==key and r['time']==1000)} for key in keys},absolute_tracer_full_field_extrema=[tmin,tmax],notes=['Zonal-u fraction uses reference-density-weighted depth-mean velocities, not total 3D kinetic energy.','Mean-flow amplitudes and stress hypotheses do not close a momentum budget or establish causation.','Both nojet cases have no imposed momentum restoration.','Tracer variance decay includes stirring and numerical diffusion; this is not a calibrated physical diffusivity.'])
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(out/'profiles.json').write_text(json.dumps(profiles)+'\n')
limits=old['vorticity_limits'];delta=max(float(np.max(abs(a-a.mean(axis=(2,3),keepdims=True,dtype=np.float64)))) for a in C)
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':8})
def write_movie(kind):
 fig,axs=plt.subplots(4,3,figsize=(15,16),dpi=90);fig.subplots_adjust(left=.055,right=.935,bottom=.045,top=.92,wspace=.60,hspace=.65)
 ims=[];texts=[]
 for i in range(4):
  for j in range(3):
   if kind=='vorticity':norm=AsinhNorm(linear_width=.01,vmin=-limits[j],vmax=limits[j]);cmap='RdBu_r';clabel='ζ/f₀';ticks=[-limits[j],-.1,0,.1,limits[j]];ticklabels=[f'{-limits[j]:.2g}','−0.1','0','0.1',f'{limits[j]:.2g}']
   elif kind=='tracer_concentration':norm=FuncNorm((lambda x:np.arcsinh((x-.1)/1e-6),lambda x:.1+1e-6*np.sinh(x)),vmin=tmin,vmax=tmax);cmap='BrBG';clabel='Concentration C';ticks=[.05,.099,.09999,.1,.10001,.101,.15];ticklabels=['.05','.099','.09999','.10000','.10001','.101','.15']
   else:norm=AsinhNorm(linear_width=1e-7,vmin=-delta,vmax=delta);cmap='BrBG';clabel='C − plane mean';ticks=[-delta,-1e-3,-1e-5,0,1e-5,1e-3,delta];ticklabels=[f'{-delta:.2g}','−10⁻³','−10⁻⁵','0','10⁻⁵','10⁻³',f'{delta:.2g}']
   im=axs[i,j].imshow(np.zeros((256,128)),origin='lower',extent=[0,4,-4,4],aspect='auto',interpolation='nearest',cmap=cmap,norm=norm);ims.append(im)
   cb=fig.colorbar(im,ax=axs[i,j],fraction=.045,pad=.035);cb.set_ticks(ticks,labels=ticklabels);cb.set_label(clabel)
   axs[i,j].set(xlabel='x/H',ylabel='y/H');texts.append(axs[i,j].set_title(''))
 title=fig.suptitle('',y=.982)
 names={'vorticity':'Vorticity ζ/f₀: positive cyclonic, negative anticyclonic','tracer_concentration':'Absolute tracer concentration · fixed asinh contrast centered at C=0.1','tracer_anomaly':'Tracer contrast C − horizontal-plane mean · fixed shared scale; mean and RMS reported'}
 def update(n):
  t=times[n]
  for i in range(4):
   for j in range(3):
    a=Q[i][n,j] if kind=='vorticity' else C[i][n,j];mu=float(a.mean(dtype=np.float64));sigma=float(a.std(dtype=np.float64))
    ims[3*i+j].set_data(a-mu if kind=='tracer_anomaly' else a)
    caption=f'{labels[i]} · z/H={heights[j]:.5f}'
    if kind!='vorticity':caption+=f'\nmean C={mu:.7f}; σC={sigma:.2e}'
    texts[3*i+j].set_text(caption)
  title.set_text(f'{names[kind]}\nt = {t} · Ω=3, 67.5°N · 128×256×80 · WENO9 / no closure · Δt_movie=5')
 for n in [0,len(times)-1]:
  update(n);fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height()
  for ax in fig.axes:
   for txt in [ax.title,ax.xaxis.label,ax.yaxis.label,*ax.get_xticklabels(),*ax.get_yticklabels()]:
    if txt.get_visible() and txt.get_text():
     b=txt.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=w and b.y1<=h,(txt.get_text(),b.bounds)
  assert w%2==h%2==0
 fig.savefig(out/f'{kind}_2000.png',dpi=90)
 movie=out/f'{kind}_0_2000.mp4';writer=FFMpegWriter(fps=8,codec='libx264',bitrate=1500,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
 with writer.saving(fig,str(movie),dpi=90):
  for n in range(len(times)):update(n);writer.grab_frame()
 plt.close(fig)
 nf,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert nf==401
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
 print(kind,'401 FRAMES FULL DECODE VERIFIED',flush=True)
for kind in ['vorticity','tracer_concentration','tracer_anomaly']:write_movie(kind)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
media=dict(status='MEDIA_VERIFIED',full_decode_passed=True,frames_per_movie=401,fps=8,times=[0,2000],sampling_interval=5,row_order=keys,heights=heights.tolist(),vorticity_limits=limits,absolute_tracer=dict(center=.1,width=1e-6,limits=[tmin,tmax]),tracer_anomaly=dict(definition='C minus instantaneous horizontal mean of the displayed plane',width=1e-7,limits=[-delta,delta],per_frame_normalization=False),assets=assets)
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)

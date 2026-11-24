from pathlib import Path
import json,csv,hashlib,subprocess,shutil,tomllib
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter

root=Path.cwd();src=root/'data/chan-clean-pair-v1-run-2214011_1/parent';out=root/'rendered/chan-clean-anelastic2000-v2';out.mkdir(parents=True,exist_ok=False)
meta=tomllib.loads((src/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==2000
times=np.arange(0,2001,5);Q=np.empty((401,3,256,256),np.float32);C=np.empty_like(Q);rows=[];lastmaps={}
ky=np.fft.fftfreq(256)*256;kx=ky;rad=np.hypot(ky[:,None],kx[None,:]);dx=4/256;f=meta['f']
for n,t in enumerate(times):
 with h5py.File(src/f'frame_{int(t*1000):09d}.jld2') as h:
  assert abs(float(h['t'][()])-t)<1e-7
  Q[n]=h['vorticity_slices'][:]/f;C[n]=h['tracer_slices'][:];heights=h['video_heights'][:]
  u=h['ubar'][:].astype(float);v=h['vbar'][:].astype(float);u-=u.mean();v-=v.mean()
  q=((np.roll(v,-1,1)-np.roll(v,1,1))-(np.roll(u,-1,0)-np.roll(u,1,0)))/(2*dx*f)
  uh=np.fft.fft2(u);vh=np.fft.fft2(v);qh=np.fft.fft2(q);spec=abs(uh)**2+abs(vh)**2;den=spec.sum()
  ke=.5*np.mean(u*u+v*v);zonal=.5*np.mean(u.mean(axis=1)**2)
  r={'time':int(t),'depthmean_nonuniform_ke':float(ke),'zonal_u_ke_fraction':float(zonal/ke) if ke>0 else 0.,'depthmean_mean_zeta_over_f':float(q.mean())}
  for k in (2,4):
   mask=(rad<=k)&(rad>0);qlo=np.fft.ifft2(qh*mask).real;en=np.sum(qlo*qlo);rms=np.sqrt(np.mean(qlo*qlo))
   r.update({f'lowmode_ke_fraction_k{k}':float(spec[mask].sum()/den) if den>0 else 0.,f'anticyclonic_enstrophy_share_k{k}':float(np.sum(qlo[qlo<0]**2)/en) if en>0 else .5,
     f'vorticity_skewness_k{k}':float(np.mean(qlo**3)/rms**3) if rms>0 else 0.,f'vorticity_min_over_f_k{k}':float(qlo.min()),f'vorticity_max_over_f_k{k}':float(qlo.max())})
   if t==2000:lastmaps[f'q{k}']=qlo
  for j in range(3):r[f'tracer_std_z{j}']=float(C[n,j].std(dtype=np.float64));r[f'tracer_mean_z{j}']=float(C[n,j].mean(dtype=np.float64))
  rows.append(r)
 if n%100==0:print('CACHED',t,flush=True)
assert np.isfinite(Q).all() and np.isfinite(C).all()
with (out/'structure_diagnostics.csv').open('w') as f0:
 w=csv.DictWriter(f0,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with (src/'diagnostics.csv').open() as f0: native=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f0)]
native_last=native[-1];summary={'case':'clean anelastic Omega=3','source_job':'2214011_1','time_end':2000,'grid':[256,256,80],
 'initial':rows[0],'final':rows[-1],'selected_times':[r for r in rows if r['time'] in (500,1000,1500,2000)],
 'native_final':native_last,'native_max_mach':max(r['max_mach'] for r in native),'native_max_constraint':max(r['constraint'] for r in native),
 'late_1500_2000':{k:float(np.mean([r[k] for r in rows if r['time']>=1500])) for k in rows[0] if k!='time'},
 'notes':['Times are nondimensional, not Earth seconds or days.','Fourier radius k is cycles per box; low-mode metrics use reference-density-weighted depth-mean horizontal velocity with its spatial mean removed.',
 'Anticyclonic share is the fraction of squared filtered vorticity in negative-vorticity cells; it is not vortex count, circulation fraction, or proof of an LSV.',
 'Periodic net vertical vorticity is zero; concentrated vorticity must have compensating circulation.','Spectral energy concentration and morphology do not prove inverse energy transfer.',
 'The compressible energy proxy is not a conserved anelastic energy invariant; its residual is not reported as numerical energy loss.',
 'The compressible partner has not yet completed, so this entry does not establish a compressibility-induced sign difference.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
shutil.copy2(src/'metadata.toml',out/'metadata.toml')
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':9,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for k,style in ((2,'-'),(4,'--')):
 axs[0].plot(times,[r[f'lowmode_ke_fraction_k{k}'] for r in rows],style,label=f'k ≤ {k}')
 axs[1].plot(times,[r[f'anticyclonic_enstrophy_share_k{k}'] for r in rows],style,label=f'k ≤ {k}')
axs[1].axhline(.5,color='gray',lw=.8);q=lastmaps['q2'];lim=np.max(abs(q));im=axs[2].imshow(q,origin='lower',extent=[0,4,0,4],cmap='RdBu_r',vmin=-lim,vmax=lim);fig.colorbar(im,ax=axs[2],label='ζ / f₀')
axs[0].set(xlabel='t',ylabel='Fraction of depth-mean nonuniform KE',title='Energy at large horizontal scales');axs[1].set(xlabel='t',ylabel='Negative-vorticity enstrophy share',title='Filtered vorticity sign statistic');axs[2].set(xlabel='x/H',ylabel='y/H',title='Depth-mean vorticity, k ≤ 2, t=2000')
for ax in axs[:2]:ax.legend();ax.grid(alpha=.2)
fig.savefig(out/'structure.svg');plt.close(fig)
qmax=np.max(abs(Q),axis=(0,2,3));cmin=float(C.min());cmax=float(C.max());late_std=float(np.median(np.std(C[300:],axis=(2,3),dtype=np.float64)));width=float(np.clip(late_std/3,1e-8,.005))
means=C.mean(axis=(2,3),keepdims=True,dtype=np.float64);anomaly_limit=max(float(np.max(abs(C[n]-means[n]))) for n in range(401))
def movie(kind):
 combo=kind=='vorticity_tracer';nr=2 if combo else 1
 fig,axs=plt.subplots(nr,3,figsize=(14,8.8 if combo else 4.8),dpi=100,squeeze=False)
 fig.subplots_adjust(left=.055,right=.94,bottom=.08,top=.87 if combo else .78,wspace=.42,hspace=.44)
 ims=[];titles=[]
 for row in range(nr):
  for j in range(3):
   if combo and row==0:norm=AsinhNorm(linear_width=.02,vmin=-qmax[j],vmax=qmax[j]);cm='RdBu_r';label='ζ/f₀'
   elif combo:norm=FuncNorm((lambda x:np.arcsinh((x-.1)/width),lambda x:.1+width*np.sinh(x)),vmin=cmin,vmax=cmax);cm='BrBG';label='Absolute C'
   else:norm=AsinhNorm(linear_width=width,vmin=-anomaly_limit,vmax=anomaly_limit);cm='BrBG';label='C − plane mean'
   im=axs[row,j].imshow(np.zeros((256,256)),origin='lower',extent=[0,4,0,4],interpolation='nearest',norm=norm,cmap=cm);ims.append(im)
   cb=fig.colorbar(im,ax=axs[row,j],fraction=.045,pad=.03,label=label)
   if combo and row==1:cb.set_ticks([cmin,.09999,.1,.10001,cmax],labels=[f'{cmin:.3f}','.09999','.10000','.10001',f'{cmax:.3f}'])
   elif not combo:cb.set_ticks([-anomaly_limit,-1e-5,-1e-7,0,1e-7,1e-5,anomaly_limit],labels=[f'{-anomaly_limit:.2g}','−10⁻⁵','−10⁻⁷','0','10⁻⁷','10⁻⁵',f'{anomaly_limit:.2g}'])
   axs[row,j].set(xlabel='x/H',ylabel='y/H');titles.append(axs[row,j].set_title('',fontsize=9))
 title=fig.suptitle('',y=.98,fontsize=14)
 def update(n):
  for row in range(nr):
   for j in range(3):
    vort=combo and row==0;a=Q[n,j] if vort else C[n,j];mu=float(C[n,j].mean(dtype=np.float64));sig=float(C[n,j].std(dtype=np.float64))
    ims[3*row+j].set_data(a if combo else a-mu)
    titles[3*row+j].set_text(('Vorticity' if vort else ('Absolute tracer' if combo else 'Tracer anomaly'))+f' · z/H={heights[j]:.5f}'+('' if vort else f'\nmean C={mu:.7f}; σ={sig:.2e}'))
  title.set_text(f'Clean anelastic convection · Ω=3, β=0, no jet · t={times[n]}\n256² × 80 · fixed color scales · '+('positive ζ/f₀ is cyclonic' if combo else 'C minus each plane’s instantaneous mean; no variance normalization'))
 update(400);fig.savefig(out/f'{kind}_2000.jpg',dpi=100,pil_kwargs={'quality':88})
 fn=out/f'{kind}_0_2000.mp4';writer=FFMpegWriter(fps=12,codec='libx264',bitrate=1600,extra_args=['-threads','2','-maxrate','2000k','-bufsize','4000k','-pix_fmt','yuv420p','-movflags','+faststart'])
 with writer.saving(fig,str(fn),dpi=100):
  for n in range(401):update(n);writer.grab_frame()
 plt.close(fig);nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==401
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
 print(kind,'FULL_DECODE_PASS',fn.stat().st_size,flush=True)
for kind in ['vorticity_tracer','tracer_anomaly']:movie(kind)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
media={'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames_per_movie':401,'fps':12,'times':[0,2000],'sampling_interval':5,
 'heights':heights.tolist(),'vorticity_limits_per_depth':qmax.tolist(),'vorticity_asinh_width':.02,
 'tracer_absolute':{'limits':[cmin,cmax],'center':.1,'asinh_width':width},'tracer_anomaly':{'limits':[-anomaly_limit,anomaly_limit],'asinh_width':width,'definition':'C minus instantaneous horizontal plane mean; no per-frame variance normalization'},'assets':assets}
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)

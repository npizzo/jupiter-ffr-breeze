"""Completed matched Omega=3 comparison; read-only computation/rendering on Unity."""
import sys,json,csv,tomllib,subprocess,hashlib
from pathlib import Path
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
a=Path('data/omega3-anelastic-long-v1-run-2113514/parent');c=Path('data/omega3-long-v1-run-2096858/parent')
old=Path('data/chan-rotation-response-v1-control-2093290/parent');cor=json.loads((old.parent/'corrected_diagnostics/CORRECTION.json').read_text())
segments={'anelastic':[Path('data/omega3-anelastic-v4-run-2113277/parent'),a],'compressible':[Path('data/chan-eddies-v3-screen-2089476/case4'),Path('data/chan-eddies-v4-case4-2090026/parent'),old,c]}
cases={};limits=np.zeros(3)
def metrics(u,v):
 n=u.shape[0];dx=4/n;k=np.fft.fftfreq(n)*n;kx,ky=np.meshgrid(k,k);kr=np.hypot(kx,ky);low=(kr>0)&(kr<=2)
 e=.5*(abs(np.fft.fft2(u)/n**2)**2+abs(np.fft.fft2(v)/n**2)**2)
 z=(np.roll(v,-1,1)-np.roll(v,1,1)-np.roll(u,-1,0)+np.roll(u,1,0))/(2*dx)
 zl=np.fft.ifft2(np.fft.fft2(z)*low).real
 skew=lambda x:float(np.mean(x**3)/max(np.mean(x*x)**1.5,1e-99))
 rec=dict(barotropic_ke=float(e.sum()),gravest_ke=float(e[np.isclose(kr,1)].sum()),low_k_le2_ke=float(e[low].sum()),low_k_anticyclonic_enstrophy_fraction=float(np.sum(zl[zl<0]**2)/np.sum(zl**2)),low_k_positive_peak=float(zl.max()),low_k_negative_peak=float(-zl.min()),barotropic_skew=skew(z),low_k_skew=skew(zl))
 return rec,z,zl
for label,segs in segments.items():
 root=segs[-1];meta=tomllib.loads((root/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==1000
 frames={};rows={}
 for seg in segs:
  for p in seg.glob('frame_*.jld2'):
   if seg==old and p.name==cor['excluded_movie_frame']:continue
   with h5py.File(p) as h:t=float(h['t'][()])
   if abs(t-round(t))>1e-8:continue
   assert t not in frames,(label,t);frames[t]=p
  dp=old.parent/'corrected_diagnostics/diagnostics_from_first_poststep.csv' if seg==old else seg/'diagnostics.csv'
  for rr in csv.DictReader(dp.open()):
   rr={k:float(v) for k,v in rr.items()};rows[rr['time']]=rr
 frames=dict(sorted(frames.items()));assert np.allclose(list(frames),np.arange(1001))
 records=[];pe=[]
 for t,p in frames.items():
  with h5py.File(p) as h:
   zz=h['vorticity_slices'][()].astype(float);tr=h['tracer_slices'][()];f=float(h['f'][()]);heights=h['video_heights'][()].tolist()
   assert zz.shape==(3,256,256) and np.isfinite(zz).all() and np.isfinite(tr).all()
   if label=='anelastic':limits=np.maximum(limits,np.max(abs(zz/f),axis=(1,2)))
   rec,_,_=metrics(h['ubar'][()].astype(float),h['vbar'][()].astype(float));rec['time']=t
   for j in range(3):rec[f'skew_{j}']=float(np.mean(zz[j]**3)/max(np.mean(zz[j]**2)**1.5,1e-99))
   err=abs(rec['barotropic_ke']/rows[t]['barotropic_ke']-1);pe.append(err);assert err<1e-6
   records.append(rec)
 windows={f'{lo}-{hi}':{k:float(np.mean([r[k] for r in records if lo<=r['time']<=hi])) for k in records[0] if k!='time'} for lo,hi in [(450,500),(700,750),(950,1000)]}
 summary=dict(final_time=1000,frame_count=1001,heights=heights,f=f,final_modes=records[-1],windows=windows,max_parseval_error=max(pe),final_diagnostics=rows[max(rows)],sources=[str(p) for p in segs])
 if label=='anelastic':
  summary['max_constraint_relative']=max(r['constraint_relative'] for r in rows.values());summary['max_mach']=max(r['max_mach'] for r in rows.values());summary['mapping']=tomllib.loads((segs[0]/'INITIAL_MAPPING.toml').read_text())
 with (out/f'{label}_modes.csv').open('w') as h:w=csv.DictWriter(h,fieldnames=records[0]);w.writeheader();w.writerows(records)
 cases[label]=dict(frames=frames,summary=summary,records=records)
 print(label,'ANALYSIS_PASS',json.dumps(windows['950-1000']),flush=True)
# Common fixed-density weighting checks at the three final full-field times.
common={};final_fields={}
for t in [960,980,1000]:
 with h5py.File(a/f'fields3d/state_{1000*t:09d}.jld2') as h:rho0=h['density'][()].astype(float)
 weight=rho0/rho0.sum(axis=0,keepdims=True)
 for label,root,name in [('anelastic',a,f'state_{1000*t:09d}.jld2'),('compressible',c,f'state_{t:07d}.jld2')]:
  with h5py.File(root/'fields3d'/name) as h:
   uf=h['u_face'][()].astype(float);vf=h['v_face'][()].astype(float)
   u=(weight*.5*(uf+np.roll(uf,-1,2))).sum(axis=0);v=(weight*.5*(vf+np.roll(vf,-1,1))).sum(axis=0)
  rec,z,zl=metrics(u,v);common[f'{label}_{t}']=rec
  if t==1000:final_fields[label]=(z,zl)
fig,axs=plt.subplots(2,2,figsize=(10,8),layout='constrained');lim=max(np.max(np.abs(v)) for pair in final_fields.values() for v in pair)
for row,label in enumerate(['compressible','anelastic']):
 for col,(field,name) in enumerate(zip(final_fields[label],['Full depth-mean ζ','Filtered 0<|k|≤2 ζ'])):
  im=axs[row,col].imshow(field,origin='lower',extent=[0,4,0,4],cmap='RdBu_r',norm=AsinhNorm(linear_width=.05,vmin=-lim,vmax=lim));axs[row,col].set_title(label.title()+' · '+name);axs[row,col].set_xlabel('x/H');axs[row,col].set_ylabel('y/H');fig.colorbar(im,ax=axs[row,col],label='ζ (cyclonic +)')
fig.suptitle('Ω=3, t=1000 · common fixed-density weighting and color scale');fig.savefig(out/'final_comparison.png',dpi=125);plt.close(fig)
fig,axs=plt.subplots(3,1,figsize=(10,9),sharex=True,layout='constrained')
for label,obj in cases.items():
 rr=obj['records'];tt=[r['time'] for r in rr]
 for ax,key in zip(axs,['low_k_anticyclonic_enstrophy_fraction','low_k_le2_ke','low_k_skew']):ax.plot(tt,[r[key] for r in rr],label=label)
for ax in axs:ax.grid(alpha=.2);ax.legend()
axs[0].axhline(.5,color='gray',ls=':');axs[0].set_ylabel('Anticyclonic share of\nfiltered enstrophy');axs[1].set_ylabel('Depth-mean low-k KE');axs[2].set_ylabel('Filtered ζ skewness');axs[2].set_xlabel('Nondimensional time')
fig.suptitle('Ω=3 · same filter (0<|k|≤2); native density weighting per model');fig.savefig(out/'evolution.png',dpi=125);plt.close(fig)
summary=dict(cases={k:v['summary'] for k,v in cases.items()},common_fixed_density_weighting=common,limitations=['Same velocity and entropy but thermodynamic pressure/density mapping differs; not an acoustics-only intervention.','One realization per equation set; no ensemble or resolution-robust causal attribution.','Anelastic dynamics retains stratification-driven expansion.','Native depth means have different density weights; common fixed-reference weighting checked at t=960,980,1000.','Filtered signed enstrophy depends on cutoff and does not uniquely identify a persistent material vortex.','Spectral occupancy is not an inverse-transfer budget.','Original compressible energy-budget residual is unresolved; no cross-model energy-closure claim.'])
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# All 1001 integer-time frames, 20 fps; fixed scales and concentration ticks.
plt.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
fig,axs=plt.subplots(2,3,figsize=(12,7));fig.subplots_adjust(left=.055,right=.94,bottom=.08,top=.88,wspace=.32,hspace=.38);ims=[]
s=cases['anelastic']['summary'];limits=np.maximum(limits,.05)
for j in range(3):
 im=axs[0,j].imshow(np.zeros((256,256)),origin='lower',extent=(0,4,0,4),cmap='RdBu_r',norm=AsinhNorm(linear_width=.02,vmin=-limits[j],vmax=limits[j]),interpolation='nearest');ims.append(im);axs[0,j].set_title(f'z/H = {heights[j]:.5f}');fig.colorbar(im,ax=axs[0,j],fraction=.046,pad=.025,label='ζ/f (cyclonic +)')
for j in range(3):
 im=axs[1,j].imshow(np.zeros((256,256)),origin='lower',extent=(0,4,0,4),cmap='BrBG',norm=AsinhNorm(linear_width=.01,vmin=0,vmax=1),interpolation='nearest');ims.append(im);cb=fig.colorbar(im,ax=axs[1,j],fraction=.046,pad=.025,label='Tracer concentration');cb.set_ticks([0,.01,.05,.1,.2,.5,1]);cb.set_ticklabels(['0','.01','.05','.1','.2','.5','1'])
for ax in axs.flat:ax.set_xlabel('x/H');ax.set_ylabel('y/H')
title=fig.suptitle('');movie=out/'anelastic_0_1000.mp4'
writer=FFMpegWriter(fps=20,codec='libx264',bitrate=700,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=90):
 for t,p in cases['anelastic']['frames'].items():
  with h5py.File(p) as h:zz=h['vorticity_slices'][()]/s['f'];tr=h['tracer_slices'][()]
  for j in range(3):ims[j].set_data(zz[j]);ims[j+3].set_data(tr[j])
  title.set_text(f'Anelastic convection · fixed Ω=3 · t={t:.0f}\nThree depths · fixed scales · nondimensional time')
  writer.grab_frame()
fig.savefig(out/'anelastic_0_1000.png',dpi=110);plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==1001 and abs(duration-50.05)<.1
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
media=dict(status='PASS',frame_count=count,duration=duration,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),bytes=movie.stat().st_size,vorticity_limits_per_depth=limits.tolist(),vorticity_asinh_width=.02,tracer_asinh_width=.01,tracer_range=[0,1],actual_heights=heights,all_integer_time_frames=True,location='Unity')
(out/'media_full.json').write_text(json.dumps(media,indent=2)+'\n');print('PUBLICATION_READY',json.dumps(media),flush=True)

from pathlib import Path
import json,csv,hashlib,subprocess,shutil,tomllib
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
from matplotlib.animation import FFMpegWriter

import sys
case=sys.argv[1]
cases={
 'anelastic10000': {'title':'Anelastic convection', 'end':10000,'sources':[(2000,'chan-clean-pair-v1-run-2214011_1'),(4000,'chan-clean-anelastic4000-v1-run-2231823'),(6000,'chan-clean-anelastic6000-v1-run-2249776'),(10000,'chan-clean-anelastic10000-v2-run-2264415')]},
 'energy-seeded2000': {'title':'Total-energy compressible · seeded eddies','end':2000,'sources':[(2000,'chan-total-energy-eddies-v1-run-2263614')]},
 'energy-rest4000': {'title':'Total-energy compressible · rest start','end':4000,'sources':[(2000,'chan-total-energy-long-v1-run-2237362'),(2050,'chan-total-energy4000-v2-run-2261539'),(4000,'chan-total-energy-resume4000-v3-run-2264203')]}}
cfg=cases[case];end=cfg['end'];root=Path.cwd();parts=[(t,root/'data'/name/'parent') for t,name in cfg['sources']];src=parts[-1][1];src0=parts[0][1];out=root/'rendered'/('chan-nov26-'+case);out.mkdir(parents=True,exist_ok=False)
meta=tomllib.loads((src/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['final_time']==end
times=np.arange(0,end+1,5);Q=np.empty((len(times),3,256,256),np.float32);C=np.empty_like(Q);rows=[];lastmaps={}
ky=np.fft.fftfreq(256)*256;kx=ky;rad=np.hypot(ky[:,None],kx[None,:]);dx=4/256;f=meta['f']
for n,t in enumerate(times):
 with h5py.File(next(part for limit,part in parts if t<=limit)/f'frame_{int(t*1000):09d}.jld2') as h:
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
   if t==end:lastmaps[f'q{k}']=qlo
  for j in range(3):r[f'tracer_std_z{j}']=float(C[n,j].std(dtype=np.float64));r[f'tracer_mean_z{j}']=float(C[n,j].mean(dtype=np.float64))
  rows.append(r)
 if n%100==0:print('CACHED',t,flush=True)
assert np.isfinite(Q).all() and np.isfinite(C).all()
with (out/'structure_diagnostics.csv').open('w') as f0:
 w=csv.DictWriter(f0,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
native=[]
for limit,part in parts:
 with (part/'diagnostics.csv').open() as f0:native.extend({k:float(v) for k,v in r.items()} for r in csv.DictReader(f0))
native_last=native[-1];summary={'case':cfg['title']+' Omega=3','source_directories':[str(p.relative_to(root)) for _,p in parts],'time_end':end,'grid':[256,256,80],
 'initial':rows[0],'final':rows[-1],'selected_times':[r for r in rows if r['time'] in (1000,2000,3000,4000,5000,6000,8000,10000)],
 'native_final':native_last,'native_max_mach':max(r['max_mach'] for r in native),'native_max_constraint':max(r['constraint'] for r in native),
 'late_final_1000':{k:float(np.mean([r[k] for r in rows if r['time']>=end-1000])) for k in rows[0] if k!='time'},
 'notes':['Times are nondimensional, not Earth seconds or days.','Fourier radius k is cycles per box; low-mode metrics use reference-density-weighted depth-mean horizontal velocity with its spatial mean removed.',
 'Anticyclonic share is the fraction of squared filtered vorticity in negative-vorticity cells; it is not vortex count, circulation fraction, or proof of an LSV.',
 'Periodic net vertical vorticity is zero; concentrated vorticity must have compensating circulation.','Spectral energy concentration and morphology do not prove inverse energy transfer.',
 'The compressible-energy proxy is not a conserved anelastic invariant; its residual is not interpreted as numerical energy loss.',
 'Sources and initial-condition differences are explicitly recorded; seeded and rest-start cases are not matched initial states.',
 'Thermal equilibrium is not established; cap stability and profile evolution are reported quantitatively.']}
# Final three-dimensional kinetic-energy distribution over horizontal wavenumbers.
# Use native staggered velocities with sqrt(face density), so Parseval recovers
# the same staggered kinetic-energy sum as the simulation diagnostic.
with h5py.File(src/f'full_{end*1000:09d}.jld2') as h:
 rho=h['density'][:].astype(np.float64)
 sp=np.zeros((256,256),float);normalization=256**4
 for field,axis in [('u_face',2),('v_face',1)]:
  vel=h[field][:].astype(np.float64);rf=.5*(rho+np.roll(rho,1,axis=axis))
  fh=np.fft.fft2(np.sqrt(rf)*vel,axes=(-2,-1));sp+=np.sum(abs(fh)**2,axis=0)/normalization
  del vel,rf,fh
 vel=h['w_face'][1:-1].astype(np.float64);rf=.5*(rho[:-1]+rho[1:])
 fh=np.fft.fft2(np.sqrt(rf)*vel,axes=(-2,-1));sp+=np.sum(abs(fh)**2,axis=0)/normalization
 total=sp.sum();scalar_ke=.5*(4*4/80)*total
 summary['final_full3d_horizontal_ke_spectrum']={'definition':'FFT of sqrt(arithmetic face density) times native staggered velocity, summed over vertical levels and all three components; Parseval preserves the staggered KE diagnostic. Saved fields are Float32.',
 'native_ke_from_saved_fields':float(scalar_ke),'native_reported_ke':native_last['kinetic_energy'],
 'fraction_k0':float(sp[rad==0].sum()/total),'fraction_k_le_2':float(sp[rad<=2].sum()/total),
 'fraction_wavelength_le_8_dx':float(sp[rad>=32].sum()/total),'fraction_wavelength_le_4_dx':float(sp[rad>=64].sum()/total),
 'limitations':'Horizontal spectrum only, one endpoint; density weighting affects spectral interpretation; does not establish vertical resolution or inverse transfer.'}
 assert abs(scalar_ke/native_last['kinetic_energy']-1)<1e-5
 shell=np.floor(rad+.5).astype(int)
 with (out/'final_horizontal_ke_spectrum.csv').open('w') as f0:
  w=csv.writer(f0);w.writerow(['horizontal_wavenumber_cycles_per_box','ke_fraction'])
  for k in range(shell.max()+1):w.writerow([k,float(sp[shell==k].sum()/total)])
 del rho,vel,rf,fh,sp
with h5py.File(src0/'frame_000000000.jld2') as h:tp0=h['mean_temperature'][:]
with h5py.File(src/f'frame_{end*1000:09d}.jld2') as h:tp1=h['mean_temperature'][:]
summary['mean_temperature_profile_departure_from_t0']=float(np.max(abs(tp1/tp0-1)))
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
shutil.copy2(src/'metadata.toml',out/'metadata.toml')
plt.rcParams.update({'animation.ffmpeg_path':imageio_ffmpeg.get_ffmpeg_exe(),'font.size':10})
qmax=np.max(abs(Q),axis=(0,2,3));cmin=float(C.min());cmax=float(C.max())
means=C.mean(axis=(-2,-1),dtype=np.float64);rms=C.std(axis=(-2,-1),dtype=np.float64)
width=float(np.clip(np.median(rms[times>=end-1000])/3,1e-8,.005))
for kind in ['vorticity','tracer']:
 isq=kind=='vorticity'
 fig,axs=plt.subplots(1,3,figsize=(16,5.3),dpi=100)
 fig.subplots_adjust(left=.05,right=.95,bottom=.13,top=.77,wspace=.52)
 images=[];titles=[]
 for j,ax in enumerate(axs):
  norm=AsinhNorm(linear_width=.02,vmin=-qmax[j],vmax=qmax[j]) if isq else FuncNorm((lambda x:np.arcsinh((x-.1)/width),lambda x:.1+width*np.sinh(x)),vmin=cmin,vmax=cmax)
  im=ax.imshow(np.zeros((256,256)),origin='lower',extent=[0,4,0,4],interpolation='nearest',cmap='RdBu_r' if isq else 'BrBG',norm=norm)
  images.append(im);titles.append(ax.set_title(''));ax.set(xlabel='x/H',ylabel='y/H')
  cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.03,label='ζ/f₀' if isq else 'Absolute C')
  if not isq:
   ticks=.1+width*np.sinh(np.linspace(np.arcsinh((cmin-.1)/width),np.arcsinh((cmax-.1)/width),5));cb.set_ticks(ticks,labels=[f'{x:.7g}' for x in ticks])
 title=fig.suptitle('',y=.97,fontsize=14)
 def update(i):
  title.set_text(f"{cfg['title']} · Ω=3, β=0, no jet · t={times[i]}\n256² × 80 · fixed color scales · "+('positive ζ/f₀ is cyclonic' if isq else 'absolute concentration; no framewise normalization'))
  for j in range(3):
   images[j].set_data(Q[i,j] if isq else C[i,j])
   titles[j].set_text(f'z/H={heights[j]:.5f}'+('' if isq else f'\nmean C={means[i,j]:.8f}; σ={rms[i,j]:.2e}'))
 update(len(times)-1);fig.savefig(out/f'{kind}_{end}.jpg',dpi=100,pil_kwargs={'quality':90})
 writer=FFMpegWriter(fps=20,codec='libx264',bitrate=1800,extra_args=['-threads','2','-maxrate','2200k','-bufsize','4400k','-pix_fmt','yuv420p','-movflags','+faststart'])
 fn=out/f'{kind}_0_{end}.mp4'
 with writer.saving(fig,str(fn),dpi=100):
  for i in range(len(times)):
   update(i);writer.grab_frame()
   if i%200==0:print(kind,'FRAME',i,flush=True)
 plt.close(fig)
 nf,duration=imageio_ffmpeg.count_frames_and_secs(str(fn));assert nf==len(times)
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(fn),'-f','null','-'],check=True)
 print(kind,'FULL_DECODE_PASS',flush=True)
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'media.json').write_text(json.dumps({'status':'MEDIA_VERIFIED','full_decode_passed':True,'frames_per_movie':len(times),'fps':20,'times':[0,end],'sampling_interval':5,'heights':heights.tolist(),'vorticity_limits_per_depth':qmax.tolist(),'vorticity_asinh_width':.02,'tracer_absolute':{'limits':[cmin,cmax],'center':.1,'asinh_width':width},'assets':assets},indent=2)+'\n')
print('ALL_MEDIA_VERIFIED',flush=True)

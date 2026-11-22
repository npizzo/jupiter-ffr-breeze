"""Read-only analysis and fixed-scale movies of seeded rotation runs. Unity only."""
import os,sys,json,csv,tomllib,hashlib,subprocess
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
base=Path('data/chan-eddies-v3-screen-2089476')
cases={};vmax=np.zeros(3);cmin=1e99;cmax=-1e99
for label,job in [('case3',2090025),('case4',2090026)]:
 root=Path(f'data/chan-eddies-v4-{label}-{job}/parent');meta=tomllib.loads((root/'metadata.toml').read_text());frames={}
 for seg in [base/label,root]:
  for p in seg.glob('frame_*.jld2'):
   with h5py.File(p) as h:t=float(h['t'][()])
   assert t not in frames,(label,t)
   frames[t]=p
 # Include all integer times, avoiding irregular 0.5 intervals of the short screen.
 frames={t:p for t,p in sorted(frames.items()) if abs(t-round(t))<1e-8}
 times=np.array(list(frames));assert times[0]==0 and np.allclose(np.diff(times),1)
 rows={}
 for seg in [base/label,root]:
  for r in csv.DictReader((seg/'diagnostics.csv').open()):
   r={k:float(v) for k,v in r.items()};rows[r['time']]=r
 records=[];pe=[]
 for t,p in frames.items():
  with h5py.File(p) as h:
   zeta=h['vorticity_slices'][()].astype(float);c=h['tracer_slices'][()];f=float(h['f'][()]);heights=h['video_heights'][()].tolist()
   assert zeta.shape==(3,256,256) and np.isfinite(zeta).all() and np.isfinite(c).all()
   vmax=np.maximum(vmax,np.max(abs(zeta/f),axis=(1,2)));cmin=min(cmin,float(c.min()));cmax=max(cmax,float(c.max()))
   u=h['ubar'][()].astype(float);v=h['vbar'][()].astype(float);n=u.shape[0];dx=4/n;k=np.fft.fftfreq(n)*n;kx,ky=np.meshgrid(k,k);kr=np.hypot(kx,ky)
   e=.5*(abs(np.fft.fft2(u)/n**2)**2+abs(np.fft.fft2(v)/n**2)**2)
   zb=(np.roll(v,-1,1)-np.roll(v,1,1)-np.roll(u,-1,0)+np.roll(u,1,0))/(2*dx)
   low=(kr>0)&(kr<=2);zl=np.fft.ifft2(np.fft.fft2(zb)*low).real
   skew=lambda a:float(np.mean(a**3)/max(np.mean(a*a)**1.5,1e-99))
   pe.append(abs(float(e.sum())-rows[t]['barotropic_ke'])/max(rows[t]['barotropic_ke'],1e-99));assert pe[-1]<1e-6
   rec=dict(time=t,barotropic_ke=float(e.sum()),gravest_ke=float(e[np.isclose(kr,1)].sum()),low_k_le2_ke=float(e[low].sum()),barotropic_skew=skew(zb),low_k_skew=skew(zl))
   for j in range(3):rec[f'skew_{j}']=skew(zeta[j])
   records.append(rec)
 ds=[rows[t] for t in sorted(rows)];last=ds[-1];tail=[r for r in records if r['time']>=times[-1]-50]
 ck=sorted((root/'checkpoints').glob('*.jld2'));assert len(ck)==1 and ck[0].stat().st_size>0
 heat=sum((b['time']-a['time'])*.5*(a['external_power']+b['external_power']) for a,b in zip(ds,ds[1:]));defect=last['total_energy']-ds[0]['total_energy']-heat
 summary=dict(job=job,omega=meta['omega'],f=f,status=meta['status'],final_time=meta['final_time'],last_saved_movie_time=float(times[-1]),frame_count=len(times),heights=heights,initial_turnovers=meta['final_time']/meta['initialization']['turnover_time'],final_diagnostics=last,initial_modes=records[0],final_modes=records[-1],last_50_time_units_means={k:float(np.mean([r[k] for r in tail])) for k in records[0] if k!='time'},energy_defect_t0_to_last_diagnostic=defect,net_heat_t0_to_last_diagnostic=heat,defect_over_net_heat=defect/heat,defect_over_initial_ke=defect/meta['initialization']['kinetic_energy'],max_mass_relative_error=max(abs(r['mass']/ds[0]['mass']-1) for r in ds),max_parseval_error=max(pe),checkpoint=str(ck[0]),checkpoint_bytes=ck[0].stat().st_size,metadata=meta)
 with (out/f'{label}_modes.csv').open('w') as h:w=csv.DictWriter(h,fieldnames=records[0]);w.writeheader();w.writerows(records)
 cases[label]=dict(frames=frames,summary=summary,records=records)
 print(label,'ANALYSIS_PASS',json.dumps({k:summary[k] for k in ['final_time','frame_count','defect_over_net_heat','final_modes','last_50_time_units_means']}),flush=True)
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ffmpeg
for label,obj in cases.items():
 fig,axs=plt.subplots(2,3,figsize=(12,7));fig.subplots_adjust(left=.055,right=.94,bottom=.08,top=.88,wspace=.32,hspace=.38)
 ims=[];s=obj['summary'];center=.1;width=.01
 for j in range(3):
  im=axs[0,j].imshow(np.zeros((256,256)),origin='lower',extent=(0,4,0,4),cmap='RdBu_r',norm=AsinhNorm(linear_width=.02,vmin=-vmax[j],vmax=vmax[j]),interpolation='nearest');ims.append(im)
  axs[0,j].set_title(f'z/H = {s["heights"][j]:.5f}');fig.colorbar(im,ax=axs[0,j],fraction=.046,pad=.025,label='ζ/f (cyclonic +)')
 for j in range(3):
  im=axs[1,j].imshow(np.zeros((256,256)),origin='lower',extent=(0,4,0,4),cmap='BrBG',norm=AsinhNorm(linear_width=width,vmin=cmin-center,vmax=cmax-center),interpolation='nearest');ims.append(im)
  cb=fig.colorbar(im,ax=axs[1,j],fraction=.046,pad=.025,label='Tracer concentration');ticks=np.linspace(cmin,cmax,5);cb.set_ticks(ticks-center);cb.set_ticklabels([f'{x:.3f}' for x in ticks])
 for ax in axs.flat:ax.set_xlabel('x/H');ax.set_ylabel('y/H')
 title=fig.suptitle('');movie=out/f'{label}.mp4'
 writer=FFMpegWriter(fps=10,codec='libx264',bitrate=1300,extra_args=['-threads','2','-pix_fmt','yuv420p','-movflags','+faststart'])
 with writer.saving(fig,str(movie),dpi=90):
  for t,p in obj['frames'].items():
   with h5py.File(p) as h:zeta=h['vorticity_slices'][()]/s['f'];c=h['tracer_slices'][()]
   for j in range(3):ims[j].set_data(zeta[j]);ims[j+3].set_data(c[j]-center)
   title.set_text(f'Seeded compressible convection · Ω = {s["omega"]:g} · t = {t:.0f}\nThree depths · fixed scales shared between cases · dimensionless time')
   writer.grab_frame()
 fig.savefig(out/f'{label}.png',dpi=110);plt.close(fig)
 subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-f','null','-'],check=True)
 s['movie']={'file':movie.name,'fps':10,'duration_seconds':s['frame_count']/10,'full_decode':'PASS','bytes':movie.stat().st_size,'sha256':hashlib.sha256(movie.read_bytes()).hexdigest()}
 print(label,'MOVIE_DECODE_PASS',movie.stat().st_size,flush=True)
summary={'analysis_job':os.getenv('SLURM_JOB_ID'),'cases':{k:v['summary'] for k,v in cases.items()},'scales':{'vorticity_limits_per_depth':vmax.tolist(),'vorticity_asinh_width':.02,'tracer_center':.1,'tracer_asinh_width':.01,'tracer_full_range':[cmin,cmax],'fixed_across_cases_and_time':True},'limitations':['Finite-amplitude seeded response, not spontaneous formation from weak noise.','One realization per rotation; no ensemble sign-transition claim.','Pressure and density balance depend on rotation; thermal background not equilibrated.','Energy budget not closed; numerical loss remains substantial.','Density-weighted depth-mean spectra are not an energy-transfer budget.','Movie includes every integer-time saved frame; half-time screen frames omitted to keep a constant playback cadence.','Last snapshot precedes final checkpoint by less than one time unit.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PUBLICATION_READY',flush=True)

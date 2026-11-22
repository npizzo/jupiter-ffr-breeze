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
for label,job in [('omega3',2096858),('omega45',2093564)]:
 if label=='omega3':
  root=Path('data/omega3-long-v1-run-2096858/parent')
  old=Path('data/chan-rotation-response-v1-control-2093290/parent')
  corrected=old.parent/'corrected_diagnostics'
  correction=json.loads((corrected/'CORRECTION.json').read_text())
  segments=[base/'case4',Path('data/chan-eddies-v4-case4-2090026/parent'),old,root]
 else:
  root=Path('data/chan-fixed45-v2-run-2093564/parent');old=None;correction=None;segments=[root]
 meta=tomllib.loads((root/'metadata.toml').read_text());frames={}
 assert meta['status']==('COMPLETE' if label=='omega3' else 'WALL_LIMIT')
 assert abs(meta['final_time']-(1000 if label=='omega3' else 491.50652840605903))<1e-7
 for seg in segments:
  for p in seg.glob('frame_*.jld2'):
   if seg==old and p.name==correction['excluded_movie_frame']:continue
   with h5py.File(p) as h:t=float(h['t'][()])
   assert t not in frames,(label,t)
   frames[t]=p
 # Include all integer times, avoiding irregular 0.5 intervals of the short screen.
 frames={t:p for t,p in sorted(frames.items()) if abs(t-round(t))<1e-8}
 times=np.array(list(frames));assert times[0]==0 and np.allclose(np.diff(times),1)
 rows={}
 for seg in segments:
  for r in csv.DictReader(((corrected/'diagnostics_from_first_poststep.csv') if seg==old else (seg/'diagnostics.csv')).open()):
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
   rec=dict(low_k_positive_peak=float(zl.max()),low_k_negative_peak=float(-zl.min()),low_k_anticyclonic_enstrophy_fraction=float(np.sum(zl[zl<0]**2)/max(np.sum(zl**2),1e-99)),time=t,barotropic_ke=float(e.sum()),gravest_ke=float(e[np.isclose(kr,1)].sum()),low_k_le2_ke=float(e[low].sum()),barotropic_skew=skew(zb),low_k_skew=skew(zl))
   for j in range(3):rec[f'skew_{j}']=skew(zeta[j])
   records.append(rec)
 ds=[rows[t] for t in sorted(rows)];last=ds[-1];tail=[r for r in records if r['time']>=times[-1]-50]
 ck=sorted((root/'checkpoints').glob('*.jld2'));assert len(ck)==1 and ck[0].stat().st_size>0
 budget_rows=[{k:float(v) for k,v in rr.items()} for rr in csv.DictReader((root/'diagnostics.csv').open())]
 heat=sum((b['time']-a['time'])*.5*(a['external_power']+b['external_power']) for a,b in zip(budget_rows,budget_rows[1:]));defect=budget_rows[-1]['total_energy']-budget_rows[0]['total_energy']-heat
 budget_interval=[budget_rows[0]['time'],budget_rows[-1]['time']]
 assert np.allclose(times,np.arange(1001 if label=='omega3' else 492))
 summary=dict(omega=meta['omega'],f=f,status=meta['status'],final_time=meta['final_time'],last_saved_movie_time=float(times[-1]),frame_count=len(times),heights=heights,initial_turnovers=meta['final_time']/meta['initialization']['turnover_time'],final_diagnostics=last,initial_modes=records[0],final_modes=records[-1],last_50_time_units_means={k:float(np.mean([r[k] for r in tail])) for k in records[0] if k!='time'},energy_budget_interval=budget_interval,energy_defect_corrected_interval=defect,net_heat_corrected_interval=heat,defect_over_net_heat=defect/heat,defect_over_initial_ke=defect/meta['initialization']['kinetic_energy'],max_mass_relative_error=max(abs(r['mass']/ds[0]['mass']-1) for r in ds),max_parseval_error=max(pe),checkpoint_present=True,restart_output_correction=correction,grid=[256,256,80],domain=[4,4,1],initial_rms=.136)
 summary['windows']={f'{a}-{b}':{k:float(np.mean([r[k] for r in records if a<=r['time']<=b])) for k in records[0] if k!='time'} for a,b in ([(213,263),(441,491),(450,500),(700,750),(950,1000)] if label=='omega3' else [(213,263),(441,491)])}
 summary['final_gravest_fraction']=records[-1]['gravest_ke']/records[-1]['barotropic_ke']
 summary['final_ke_over_initial']=records[-1]['barotropic_ke']/records[0]['barotropic_ke']
 summary['last50_skew_ranges']={k:[min(r[k] for r in tail),max(r[k] for r in tail)] for k in ['barotropic_skew','low_k_skew','skew_0','skew_1','skew_2']}
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
   title.set_text(f'Seeded compressible convection · Ω = {s["omega"]:g} · t = {t:.0f}\nThree depths · fixed scales shared across both runs · dimensionless time')
   writer.grab_frame()
 fig.savefig(out/f'{label}.png',dpi=110);plt.close(fig)
 subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-f','null','-'],check=True)
 s['movie']={'file':movie.name,'fps':10,'duration_seconds':s['frame_count']/10,'full_decode':'PASS','bytes':movie.stat().st_size,'sha256':hashlib.sha256(movie.read_bytes()).hexdigest()}
 print(label,'MOVIE_DECODE_PASS',movie.stat().st_size,flush=True)
summary={'analysis_location':'Unity','cases':{k:v['summary'] for k,v in cases.items()},'scales':{'vorticity_limits_per_depth':vmax.tolist(),'vorticity_asinh_width':.02,'tracer_center':.1,'tracer_asinh_width':.01,'tracer_full_range':[cmin,cmax],'fixed_across_cases_and_time':True},'limitations':['Finite-amplitude seeded response, not spontaneous formation from weak noise.','One realization per rotation; no ensemble sign-transition claim.','Pressure and density balance depend on rotation; thermal background not equilibrated.','Energy budget not closed; numerical loss remains substantial.','Density-weighted depth-mean spectra are not an energy-transfer budget.','Movie includes every integer-time saved frame; half-time screen frames omitted to keep a constant playback cadence.','Old Omega=3 invalid restart sample omitted and correction preserved. New continuation budget starts at its first refreshed output, t=500.5. Omega=4.5 budget starts at t=0.', 'Omega=4.5 stopped at the configured wall limit at t=491.5065; its last integer-time movie frame is t=491. This is partial relative to target t=500.', 'Low-pass peak and signed enstrophy diagnostics depend on the specified cutoff and do not track individual coherent vortex identities.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for label,obj in cases.items():
 rr=obj['records'];tt=[r['time'] for r in rr]
 fig,axs=plt.subplots(4,1,figsize=(10,11),sharex=True)
 for k,name in [('barotropic_ke','All depth-mean KE'),('gravest_ke','Gravest |k|=1 KE'),('low_k_le2_ke','0<|k|≤2 KE')]:axs[0].plot(tt,[r[k] for r in rr],label=name)
 for k,name in [('barotropic_skew','Full depth-mean ζ'),('low_k_skew','Filtered large-scale ζ')]:axs[1].plot(tt,[r[k] for r in rr],label=name)
 for j,z in enumerate(obj['summary']['heights']):axs[2].plot(tt,[r[f'skew_{j}'] for r in rr],label=f'z/H={z:.5f}')
 for k,name in [('low_k_positive_peak','Cyclonic peak'),('low_k_negative_peak','Anticyclonic peak magnitude')]:axs[3].plot(tt,[r[k] for r in rr],label=name)
 axs[0].set_ylabel('Kinetic energy');axs[1].set_ylabel('Vorticity skewness');axs[2].set_ylabel('Slice skewness');axs[3].set_ylabel('Filtered peak |ζ|');axs[3].set_xlabel('Nondimensional time')
 for ax in axs:
  if label=='omega3':
   for t in [263.06682088232816,500]:ax.axvline(t,color='gray',ls=':')
  ax.axhline(0,color='gray',lw=.5);ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.suptitle(f'Constant Ω={obj["summary"]["omega"]:g} · depth-mean KE and sign diagnostics')
 fig.tight_layout();fig.savefig(out/f'{label}_evolution.png',dpi=130);plt.close(fig)
 # Diagnostic illustration: fixed common signed limits within this final-frame figure.
 with h5py.File(list(obj['frames'].values())[-1]) as h:
  u=h['ubar'][()].astype(float);v=h['vbar'][()].astype(float);n=u.shape[0];dx=4/n
  zb=(np.roll(v,-1,1)-np.roll(v,1,1)-np.roll(u,-1,0)+np.roll(u,1,0))/(2*dx)
  k=np.fft.fftfreq(n)*n;kx,ky=np.meshgrid(k,k);kr=np.hypot(kx,ky)
  zl=np.fft.ifft2(np.fft.fft2(zb)*((kr>0)&(kr<=2))).real
 fig,axs=plt.subplots(1,2,figsize=(10,4));lim=max(abs(zb).max(),abs(zl).max())
 for ax,field,name in zip(axs,[zb,zl],['Full depth-mean ζ','Filtered 0<|k|≤2 ζ']):
  im=ax.imshow(field,origin='lower',extent=(0,4,0,4),cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_title(name);ax.set_xlabel('x/H');ax.set_ylabel('y/H');fig.colorbar(im,ax=ax,label='ζ')
 fig.suptitle(f'Ω={obj["summary"]["omega"]:g} · t={obj["summary"]["last_saved_movie_time"]:g} · positive = cyclonic')
 fig.tight_layout();fig.savefig(out/f'{label}_depthmean.png',dpi=130);plt.close(fig)
print('PUBLICATION_READY',flush=True)

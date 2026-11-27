from pathlib import Path
import json,shutil,hashlib,subprocess
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
root=Path.cwd();out=root/'rendered/deformation100-contrast-v2';out.mkdir(exist_ok=False,parents=True)
cases={'weak':2311446,'strong':2311447};cache={};qlim=np.zeros(3);clo=np.inf;chi=-np.inf
summary={'scope':'Unforced, stably stratified anelastic evolution of identical mixed-sign velocity fields. Internal deformation; not convection-driven LSV formation or a compressible Lamb-mode test.', 'time_note':'Nondimensional time. Integration reached t=10000 (100 initial turnovers). Original 25-turnover snapshot missing: the movie preserves the gap from 24.9 to 25.1 without interpolation. Restart-segment scalar mass errors are relative to each segment start.', 'resolution':[128,128,64],'domain':[8,8,1],'f':.1,'Ro':.1,'Re':2000,'seed':261130,'cases':{}}
for name,job in cases.items():
 p=root/f'data/deformation-continue-v1-{name}-{job}/parent';old=2308454 if name=='weak' else 2308455;parts=[root/f'data/deformation-pair-v1-{name}-{old}/parent',p.parent/'pilot',p];files=sorted([f for part in parts for f in part.glob('frame_*.jld2')],key=lambda f:f.name);files=list({f.name:f for f in files}.values());a=[];b=[];ts=[]
 for file in files:
  with h5py.File(file) as h:
   a.append(h['vorticity_slices'][:]/float(h['f'][()]));b.append(h['tracer_slices'][:]);ts.append(float(h['turnovers'][()]));heights=h['video_heights'][:]
 a=np.array(a);b=np.array(b);assert a.shape[1:]==(3,128,128) and len(ts)>=1000 and np.isfinite(a).all() and np.isfinite(b).all()
 qlim=np.maximum(qlim,abs(a).max(axis=(0,2,3)));clo=min(clo,float(b.min()));chi=max(chi,float(b.max()))
 d=np.concatenate([np.genfromtxt(part/'diagnostics.csv',delimiter=',',names=True) for part in parts]);d=d[np.argsort(d['time'])];np.savetxt(out/f'{name}_diagnostics.csv',np.column_stack([d[k] for k in d.dtype.names]),delimiter=',',header=','.join(d.dtype.names),comments='')
 stats=[]
 for i,t in enumerate(ts):
  q=a[i].astype(float);share=np.sum(np.minimum(q,0)**2,axis=(1,2))/np.sum(q*q,axis=(1,2));skew=np.mean(q**3,axis=(1,2))/np.mean(q*q,axis=(1,2))**1.5
  stats.append([t,*share,*skew])
 np.savetxt(out/f'{name}_plane_sign.csv',stats,delimiter=',',header='turnovers,anti_low,anti_mid,anti_upper,skew_low,skew_mid,skew_upper',comments='')
 q=a[-1].astype(float);rec={'job':job,'L_over_Ld':.5 if name=='weak' else 2.,'frames_saved':len(ts),'saved_last_turnovers':ts[-1],'integration_final_turnovers':100,'heights':heights.tolist(),'KE_initial':float(d['kinetic_energy'][0]),'KE_last':float(d['kinetic_energy'][-1]),'APE_proxy_initial':float(d['ape_proxy'][0]),'APE_proxy_last':float(d['ape_proxy'][-1]),'KE_plus_APE_proxy_fractional_change':float((d['kinetic_energy'][-1]+d['ape_proxy'][-1])/(d['kinetic_energy'][0]+d['ape_proxy'][0])-1),'max_mach':float(d['max_mach'].max()),'max_constraint':float(d['constraint'].max()),'min_local_N2':float(d['min_N2'].min()),'max_scalar_mass_error':float(max(d['theta_mass_relative_error'].max(),d['tracer_mass_relative_error'].max())),'upper_anti_enstrophy_last':stats[-1][3],'upper_skew_last':stats[-1][6],'last_sign_stats_all_depths':stats[-1][1:],'late_upper_anti_mean_last5turnovers':float(np.array(stats)[-50:,3].mean()),'source':str(p.relative_to(root)),'caveat':'KE + small-amplitude APE is a proxy, not an exact anelastic invariant. Changing N changes initial APE and susceptibility to baroclinic instability.'}
 summary['cases'][name]=rec;cache[name]=(a,b,np.array(ts),heights,d,np.array(stats));print(name,json.dumps(rec),flush=True)
summary['display']={'vorticity':'unfiltered cell-centered zeta/f; fixed symmetric asinh width 0.02 with full limits per depth, shared across both cases','vorticity_limits':qlim.tolist(),'tracer':'absolute concentration, fixed asinh width 0.01; full measured global range, shared across time, depths and cases','tracer_range':[clo,chi],'fps':20,'sample_interval_turnovers':.1}
from matplotlib.colors import FuncNorm
late_rms=[float(np.std(b[-100:],axis=(-2,-1)).mean()) for a,b,*_ in cache.values()]
width=max(1e-8,min(late_rms)/3)
summary['display']['tracer']='Absolute concentration, fixed centered asinh about C=0.1; width chosen from the smaller late-time RMS/3 across cases; full global extrema. Mean/RMS shown.'
summary['display']['tracer_width']=width
for name,(a,b,ts,heights,d,stats) in cache.items():
 fig,axs=plt.subplots(2,3,figsize=(14,8),dpi=100);fig.subplots_adjust(left=.06,right=.94,bottom=.08,top=.86,hspace=.4,wspace=.42);images=[]
 for j in range(3):
  im=axs[0,j].imshow(a[0,j],origin='lower',extent=[0,8,0,8],cmap='RdBu_r',norm=AsinhNorm(linear_width=.02,vmin=-qlim[j],vmax=qlim[j]),interpolation='nearest');images.append(im);fig.colorbar(im,ax=axs[0,j],fraction=.045,pad=.03,label='ζ/f (+ cyclonic)')
  im=axs[1,j].imshow(b[0,j],origin='lower',extent=[0,8,0,8],cmap='viridis',norm=FuncNorm((lambda x:np.arcsinh((x-.1)/width),lambda x:.1+width*np.sinh(x)),vmin=clo,vmax=chi),interpolation='nearest');images.append(im);cb=fig.colorbar(im,ax=axs[1,j],fraction=.045,pad=.03,label='Tracer C (absolute)');ticks=.1+width*np.sinh(np.linspace(np.arcsinh((clo-.1)/width),np.arcsinh((chi-.1)/width),5));cb.set_ticks(ticks,labels=[f'{v:.7g}' for v in ticks])
  for r in range(2):axs[r,j].set(xlabel='x/H',ylabel='y/H',title=f'z/H = {heights[j]:.4f}')
 title=fig.suptitle('',fontsize=13);fig.text(.5,.012,'128² × 64 · unforced anelastic · β = 0 · no jet · fixed colors, no temporal interpolation',ha='center',fontsize=10)
 movie=out/f'{name}_vorticity_tracer.mp4';writer=imageio_ffmpeg.write_frames(str(movie),(1400,800),fps=20,codec='libx264',bitrate='1200k',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-threads','2','-movflags','+faststart']);writer.send(None)
 for i,t in enumerate(ts):
  title.set_text(f'{name.capitalize()} deformation · L/L_D = {summary["cases"][name]["L_over_Ld"]:g} · {t:.1f} turnover times\nVorticity (top) and passive tracer (bottom)')
  for j in range(3):images[2*j].set_data(a[i,j]);images[2*j+1].set_data(b[i,j]);axs[1,j].set_title(f'z/H={heights[j]:.4f} · mean={b[i,j].mean():.6f}\nRMS variation={b[i,j].std():.2e}')
  fig.canvas.draw();writer.send(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
 writer.close();fig.savefig(out/f'{name}_final.jpg',dpi=100);plt.close(fig)
 frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==len(ts)
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
 summary['cases'][name]['media']={'frames':frames,'duration':duration,'full_decode_passed':True,'sha256':hashlib.sha256(movie.read_bytes()).hexdigest()}
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for name,(_,_,ts,_,d,s) in cache.items():
 axs[0].plot(ts,s[:,3],label=name);axs[1].plot(d['turnovers'],d['kinetic_energy'],label=name);axs[2].plot(d['turnovers'],d['ape_proxy'],label=name)
axs[0].axhline(.5,c='gray',ls='--');axs[0].set(ylabel='Upper-plane anticyclonic enstrophy share',ylim=(0,1));axs[1].set(ylabel='Volume-integrated kinetic energy');axs[2].set(ylabel='Small-amplitude APE proxy')
for ax in axs:ax.set_xlabel('Initial turnover times');ax.legend();ax.grid(alpha=.2)
fig.savefig(out/'diagnostics.png',dpi=130);plt.close(fig)
summary['status']='MEDIA_VERIFIED';(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)

from pathlib import Path
import json,tomllib,shutil,subprocess,hashlib
import numpy as np,h5py,imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm,FuncNorm
root=Path.cwd();out=root/'rendered/convective-deformation-v2';out.mkdir(parents=True,exist_ok=False)
assessment=json.loads((root/'runs/convective-deformation-v2/pilot-assessment.json').read_text())
# Preserve the original audit; correct its inherited wording only in this publication.
for r in assessment:r['note']='Cap retention measured in 0.7<z/H<0.9. V2 uses horizontal-mean cap restoration, not the old top cooling slab. Spectra at t=0 describe roundoff-level w and have no physical significance.'
summary={'status':'RENDERING','scope':'Completed short convection/stable-layer pilots. Thermal repair passes; horizontal convection-resolution screening fails. No long extension.','grid':[256,256,80],'domain':[8,8,1],'f':.1,'initial_Ro':.1,'eddy_Re':100,'Pr':1,'seed':261130,'assessment':assessment,'cases':{}}
cache={};qlim=np.zeros(3);clo=np.inf;chi=-np.inf
for label,job in [('weak',2320449),('strong',2320450)]:
 p=root/f'data/convective-deformation-v2-{label}-{job}/parent';meta=tomllib.loads((p/'metadata.toml').read_text());assert meta['status']=='COMPLETE'
 a=[];b=[];ts=[];native=[]
 for f in sorted(p.glob('frame_*.jld2')):
  with h5py.File(f) as h:a.append(h['vorticity_slices'][:]/float(h['f'][()]));b.append(h['tracer_slices'][:]);ts.append(float(h['turnovers'][()]));native.append(float(h['t'][()]));heights=h['video_heights'][:]
 a=np.array(a);b=np.array(b);ts=np.array(ts);assert a.shape==(101,3,256,256) and np.isfinite(a).all() and np.isfinite(b).all();assert np.allclose(ts,np.arange(101)/10)
 qlim=np.maximum(qlim,abs(a).max(axis=(0,2,3)));clo=min(clo,float(b.min()));chi=max(chi,float(b.max()))
 d=np.genfromtxt(p/'diagnostics.csv',delimiter=',',names=True);shutil.copy2(p/'diagnostics.csv',out/f'{label}_diagnostics.csv');shutil.copy2(p/'metadata.toml',out/f'{label}_metadata.toml')
 r=next(r for r in assessment if label in r['source']);heat=d['cumulative_bottom_heat'][-1];cap=d['cumulative_cap_heat'][-1]
 summary['cases'][label]={'job':job,'L_over_Ld':.5 if label=='weak' else 2.,'source':str(p.relative_to(root)),'frames':len(ts),'final_time':native[-1],'final_turnovers':ts[-1],'heights':heights.tolist(),'minimum_core_cap_N2_fraction':float(np.min(d['interior_cap_min_mean_N2'])/r['initial_N2']),'maximum_Mach':float(d['max_mach'].max()),'cumulative_bottom_heat':float(heat),'cumulative_cap_heat':float(cap),'net_cap_heat_over_bottom_heat':float(cap/heat),'final_cap_positive_power_over_bottom':float(d['cap_positive_power'][-1]/d['low_heat_power'][-1]),'final_cap_negative_power_over_bottom':float(d['cap_negative_power'][-1]/d['low_heat_power'][-1]),'final_w_energy_wavelength_le_8dx':r['snapshots'][-1]['w_energy_wavelength_le_8dx'],'final_peak_w_wavelength_cells':r['snapshots'][-1]['peak_w_wavelength_cells'],'production_seconds':r['wall_seconds'],'source_accounting_note':'Cap power is cp*Pi times the direct rho-theta source. Cumulative cap heat uses per-step trapezoids, not an exact total-energy closure. The cap can heat and cool at different heights; net alone can hide cancellation.'}
 cache[label]=(a,b,ts,heights,d)
width=max(1e-8,min(float(np.std(b[-20:],axis=(-2,-1)).mean()) for a,b,*_ in cache.values())/3)
summary['display']={'vorticity_asinh_width':.02,'vorticity_limits_per_depth':qlim.tolist(),'tracer_center':.1,'tracer_asinh_width':width,'tracer_limits':[clo,chi],'fps':10,'normalization':'Fixed full-extrema limits shared across both cases and all times. Absolute tracer, no mean subtraction or per-frame normalization.'}
for label,(a,b,ts,z,d) in cache.items():
 fig,axs=plt.subplots(2,3,figsize=(14,8),dpi=100);fig.subplots_adjust(left=.06,right=.94,bottom=.08,top=.86,hspace=.42,wspace=.42);images=[]
 for j in range(3):
  im=axs[0,j].imshow(a[0,j],origin='lower',extent=[0,8,0,8],cmap='RdBu_r',norm=AsinhNorm(linear_width=.02,vmin=-qlim[j],vmax=qlim[j]),interpolation='nearest');images.append(im);fig.colorbar(im,ax=axs[0,j],fraction=.045,pad=.03,label='ζ/f (+ cyclonic)')
  im=axs[1,j].imshow(b[0,j],origin='lower',extent=[0,8,0,8],cmap='BrBG',norm=FuncNorm((lambda x:np.arcsinh((x-.1)/width),lambda x:.1+width*np.sinh(x)),vmin=clo,vmax=chi),interpolation='nearest');images.append(im);cb=fig.colorbar(im,ax=axs[1,j],fraction=.045,pad=.03,label='Absolute tracer C');ticks=.1+width*np.sinh(np.linspace(np.arcsinh((clo-.1)/width),np.arcsinh((chi-.1)/width),5));cb.set_ticks(ticks,labels=[f'{v:.6g}' for v in ticks])
  for ax in axs[:,j]:ax.set(xlabel='x/H',ylabel='y/H',title=f'z/H={z[j]:.5f}')
 title=fig.suptitle('',fontsize=13);fig.text(.5,.012,'256² × 80 · dry anelastic · bottom heating + maintained mean cap · β=0 · no jet · resolution screening failed',ha='center',fontsize=10)
 movie=out/f'{label}_vorticity_tracer.mp4';writer=imageio_ffmpeg.write_frames(str(movie),(1400,800),fps=10,codec='libx264',bitrate='1800k',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-threads','2','-movflags','+faststart']);writer.send(None)
 for i,t in enumerate(ts):
  title.set_text(f'{label.capitalize()} deformation · initial L/L_D={summary["cases"][label]["L_over_Ld"]:g} · t={100*t:.0f} ({t:.1f} initial turnovers)\nVorticity above, passive tracer below · fixed scales shared between cases')
  for j in range(3):images[2*j].set_data(a[i,j]);images[2*j+1].set_data(b[i,j]);axs[1,j].set_title(f'z/H={z[j]:.5f} · mean={b[i,j].mean():.6f}\nRMS variation={b[i,j].std():.2e}')
  fig.canvas.draw();writer.send(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy())
 writer.close();fig.savefig(out/f'{label}_final.jpg',dpi=100);plt.close(fig)
 n,secs=imageio_ffmpeg.count_frames_and_secs(str(movie));assert n==101
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
 summary['cases'][label]['media']={'frames':n,'duration':secs,'full_decode_passed':True,'sha256':hashlib.sha256(movie.read_bytes()).hexdigest()};print(label,'MEDIA_PASS',flush=True)
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for label,(_,_,ts,_,d) in cache.items():
 r=next(r for r in assessment if label in r['source']);axs[0].plot(d['turnovers'],d['interior_cap_min_mean_N2']/r['initial_N2'],label=label);rr=r['snapshots'][1:];axs[1].plot([q['time']/100 for q in rr],[q['w_energy_wavelength_le_8dx'] for q in rr],'-o',label=label);axs[2].plot(d['turnovers'],d['cap_relaxation_power']/d['low_heat_power'],label=label)
axs[0].set(ylabel='Minimum mean cap N² / initial N²');axs[1].set(ylabel='Fraction of lower w power at λ ≤ 8Δx',ylim=(0,1));axs[1].axhline(.05,c='gray',ls='--');axs[2].set(ylabel='Net cap source power / bottom heat')
for ax in axs:ax.set_xlabel('Initial turnover times');ax.legend();ax.grid(alpha=.2)
fig.savefig(out/'screening.png',dpi=130);plt.close(fig)
summary['status']='MEDIA_VERIFIED';(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)

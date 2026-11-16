"""Unity-only comparison with unchanged gate thresholds; always retain failure report."""
from pathlib import Path
import sys,json,tomllib,csv,hashlib,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg
coarse,fine,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
results={};data={}
for name,p in [('256',coarse),('512',fine)]:
 c=tomllib.loads((p/'config.toml').read_text());r=tomllib.loads((p/'review.toml').read_text())
 assert r['status']=='COMPLETED'
 assert c['viscosity']==.0001 and c['cubic_velocity_drag']==0 and c['amplitude']==0
 assert c['dt']==.0025 and c['f']==1 and c['planetary_beta']==0
 a=np.genfromtxt(p/'statistics.csv',delimiter=',',names=True)
 assert np.all(np.isfinite(a.view(float)))
 s=np.genfromtxt(p/'spectra.csv',delimiter=',',names=True);nk=len(np.unique(s['k']));S={k:s[k].reshape(-1,nk) for k in s.dtype.names};k=S['k'][0]
 late=S['t'][:,0]>=4;flux=S['Pi_advection'][late]
 means=[float(flux[:,k==K].mean()) for K in [4,6,8]]
 neg=[float((flux[:,k==K]<0).mean()) for K in [4,6,8]]
 low=S['energy'][:,k<=4].sum(axis=1);gain=float((low[-1]-low[0])/a['energy'][0])
 frames=list(csv.DictReader((p/'frames.csv').open()));ts=np.array([float(x['t']) for x in frames])
 assert len(frames)==17 and np.allclose(ts,np.arange(17)*.5)
 n=int(name);kx=np.fft.fftfreq(n)*n;k2=kx[:,None]**2+kx[None,:]**2
 tails=[];fixed=[]
 def field(row):return np.fromfile(p/row['file'],dtype=np.float32).reshape((n,n),order='F').astype(float)
 for row in frames:
  q=field(row);assert np.isfinite(q).all();power=abs(np.fft.fft2(q))**2
  tails.append(float(power[k2>=(n/3)**2].sum()/power.sum()))
  fixed.append(float(power[k2>=(256/3)**2].sum()/power.sum()))
 closure=abs(a['energy'][-1]-a['energy'][0]-sum(r['energy_increments']))/a['energy'][0]
 checks={'negative_mean_flux':max(means)<0,'negative_flux_fraction_at_least_75_percent':min(neg)>=.75,'absolute_low_k_gain_at_least_one_percent_E0':gain>=.01,'high_k_enstrophy_below_one_percent':max(tails)<.01,'energy_budget_relative_error_below_1e_minus_8':closure<1e-8,'zero_forcing_and_drag_work':abs(r['energy_increments'][1])+abs(r['energy_increments'][3])<1e-20,'CFL_below_half':r['max_cfl']<.5}
 checks={key:bool(value) for key,value in checks.items()}
 results[name]={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'config':c,'review':r,'mean_flux_K4_K6_K8_t_ge_4':means,'negative_fraction_K4_K6_K8':neg,'spectral_time_endpoints':[float(S['t'][0,0]),float(S['t'][-1,0])],'first_last_spectral_low_k_energy':[float(low[0]),float(low[-1])],'low_k_gain_over_initial_energy':gain,'max_enstrophy_fraction_above_N_over_3':max(tails),'enstrophy_fraction_above_N_over_3':tails,'enstrophy_fraction_above_fixed_k_256_over_3':fixed,'relative_energy_budget_error':float(closure)}
 data[name]=(p,c,a,S,k,frames,ts)
# Refinement alone moving the tail cutoff is not evidence of convergence.
x=np.array(results['256']['mean_flux_K4_K6_K8_t_ge_4']);y=np.array(results['512']['mean_flux_K4_K6_K8_t_ge_4'])
comparison={'relative_mean_flux_changes_K4_K6_K8':((y-x)/np.maximum(abs(x),1e-30)).tolist(),'initialization':tomllib.loads((fine.parent/'initialization.toml').read_text()),'interpretation':'A tail-screen pass is not convergence. Compare fluxes, common-wavenumber tails and large-scale evolution before choosing production resolution.'}
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for name in ['256','512']:
 p,c,a,S,k,frames,ts=data[name]
 axs[0,0].plot(a['t'],a['energy'],label=name+'²')
 axs[0,1].plot(S['t'][:,0],S['energy'][:,k<=4].sum(1),label=name+'²')
 axs[1,0].semilogy(ts,np.maximum(results[name]['enstrophy_fraction_above_N_over_3'],1e-18),label=name+'²: k ≥ N/3')
 axs[1,0].semilogy(ts,np.maximum(results[name]['enstrophy_fraction_above_fixed_k_256_over_3'],1e-18),'--',label=name+'²: k ≥ 85.33')
 for K in [4,6,8]:axs[1,1].plot(S['t'][:,0],S['Pi_advection'][:,k==K],label=f'{name}², K={K}')
axs[0,0].set(title='Total energy',ylabel='E');axs[0,1].set(title='Absolute energy at |k| ≤ 4',ylabel='E_low')
axs[1,0].set(title='Enstrophy tails: moving and common cutoffs',ylabel='Fraction');axs[1,0].axhline(.01,c='k',lw=.7)
axs[1,1].set(title='Native flux: negative means upscale',ylabel='Π(K)');axs[1,1].axhline(0,c='k',lw=.7)
for ax in axs.flat:ax.set_xlabel('Model time t');ax.legend(fontsize=7)
fig.savefig(out/'resolution_comparison.png',dpi=140);plt.close(fig)
def readfield(name,row):
 p=data[name][0];n=int(name);return np.fromfile(p/row['file'],dtype=np.float32).reshape((n,n),order='F')
limit=max(float(abs(readfield(name,row)).max()) for name in data for row in data[name][5])
ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
fig,axs=plt.subplots(1,2,figsize=(12,5.5),layout='constrained');images=[]
for name,ax in zip(['256','512'],axs):
 images.append(ax.imshow(readfield(name,data[name][5][0]).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),cmap='RdBu_r',norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit),interpolation='nearest'))
 ax.set(title=name+'² · unforced · viscosity only',xlabel='x',ylabel='y')
fig.colorbar(images[0],ax=axs,label='Vertical vorticity ζ · fixed shared asinh scale',shrink=.85)
title=fig.suptitle('');fig.supxlabel('Same initial vorticity Fourier modes · dt=0.0025 · ν=10⁻⁴ · no drag\nOne horizontal 2D plane; positive ζ is counterclockwise',fontsize=10)
movie=out/'unforced_256_vs_512.mp4';writer=FFMpegWriter(fps=2,codec='libx264',extra_args=['-crf','23','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=100):
 for i,t in enumerate(data['256'][6]):
  for name,im in zip(['256','512'],images):im.set_data(readfield(name,data[name][5][i]).T)
  title.set_text(f'Unforced spatial-resolution check · t = {t:.1f} / 8');writer.grab_frame()
  if i==16:fig.savefig(out/'poster.png',dpi=100)
plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==17
subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
report={'status':results['512']['status'],'compute_location':'Unity','cases':results,'comparison':comparison,'movie':{'file':movie.name,'frames':count,'duration_seconds':duration,'full_decode_passed':True},'color_scale':{'type':'symmetric asinh','linear_width':.1,'limits':[-limit,limit],'fixed_shared':True,'clipping':False},'production_authorized_by_this_check':False,'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}}
(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':report['status'],'fine_checks':results['512']['checks'],'fine_max_tail':results['512']['max_enstrophy_fraction_above_N_over_3'],'comparison':comparison},indent=2),flush=True)
if report['status']!='PASS':sys.exit(2)

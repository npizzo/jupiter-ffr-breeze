"""Unity-only viscosity-only comparison, diagnostics and synchronized vorticity movie."""
from pathlib import Path
import csv,json,sys,tomllib,hashlib,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg
src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
cases=['forced','unforced','random_sign','cyclonic_sign']
labels=['Weak instability','Unforced decay','Weak random-sign stirring','Weak cyclonic-core stirring']
data={};report={'status':'COMPLETED_REQUIRES_REVIEW','compute_location':'Unity','cases':{},'limitations':['Finite transient comparison, not a stationary inertial-range or convergence claim.','Uniform f cannot itself select vorticity sign in this horizontal incompressible 2D model.','Stochastic core polarity is imposed; its transfer to large scales is the test.','Movies alternate post-pulse and mid-cycle states; statistics sweep all discrete pulse phases.']}
for case in cases:
 p=src/case;c=tomllib.loads((p/'config.toml').read_text());r=tomllib.loads((p/'review.toml').read_text())
 assert r['status']=='COMPLETED' and c['cubic_velocity_drag']==0 and c['viscosity']==.0001
 a=np.genfromtxt(p/'statistics.csv',delimiter=',',names=True)
 s=np.genfromtxt(p/'spectra.csv',delimiter=',',names=True);k=np.unique(s['k']);nk=len(k)
 S={key:s[key].reshape(-1,nk) for key in s.dtype.names}
 assert np.all(S['k']==k[None,:]) and np.all(np.diff(S['iteration'][:,0])>0)
 frames=list(csv.DictReader((p/'frames.csv').open()));ts=np.array([float(row['t']) for row in frames])
 assert len(frames)==401 and np.allclose(ts,np.arange(401)*.5)
 assert abs(a['energy'][-1]-a['energy'][0]-sum(r['energy_increments']))/a['energy'][0]<1e-8
 assert abs(r['energy_increments'][3])<1e-20
 late=np.flatnonzero(S['t'][:,0]>=100);late=late[:len(late)//80*80]
 assert len(late)>0
 if case in ['random_sign','cyclonic_sign']:
  phases=S['iteration'][late,0].astype(int)%80;counts=np.bincount(phases,minlength=80);assert np.all(counts==counts[0])
  e=np.genfromtxt(p/'events.csv',delimiter=',',names=True);assert len(e)==1000
  assert np.all(e['core_sign']==1) if case=='cyclonic_sign' else set(e['core_sign'])=={-1,1}
  assert np.isclose(e['energy_total'].sum(),r['energy_increments'][1],atol=1e-10)
 ia=int(np.argmin(abs(a['t']-100)));ib=len(a)-1;duration=a['t'][ib]-a['t'][ia]
 work={key:float((a['energy_'+key][ib]-a['energy_'+key][ia])/duration) for key in ['advection','injection','viscosity','drag']}
 storage=float((a['energy'][ib]-a['energy'][ia])/duration);assert abs(storage-sum(work.values()))<1e-10
 # Match phase-balanced spectral sample iterations to statistics rows.
 indices=np.searchsorted(a['iteration'],S['iteration'][late,0]);assert np.array_equal(a['iteration'][indices],S['iteration'][late,0])
 report['cases'][case]={'config':c,'review':r,'late_sampled_mean_flux_K2_to_K10':[float(S['Pi_advection'][late][:,k==K].mean()) for K in range(2,11)],'late_sample_count':len(late),'late_mean_full_skewness':float(a['skewness'][indices].mean()),'late_mean_lowpass_skewness':float(a['lowpass_skewness'][indices].mean()),'late_mean_box_energy_fraction':float(a['low_k_energy_fraction'][indices].mean()),'energy_initial_final':[float(a['energy'][0]),float(a['energy'][-1])],'late_integrated_work_rates':work,'actual_budget_window':[float(a['t'][ia]),float(a['t'][ib])],'energy_storage_rate':storage,'energy_balance_residual':storage-sum(work.values())}
 data[case]=(p,c,r,a,S,k,frames,ts,late)
assert len({data[c][1]['source_checkpoint_sha256'] for c in cases})==1
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for case,label in zip(cases,labels):
 p,c,r,a,S,k,frames,ts,late=data[case]
 axs[0,0].plot(a['t'],a['energy'],label=label)
 axs[0,1].plot(a['t'],a['energy']*a['low_k_energy_fraction'],label=label)
 for K,ls in [(4,'-'),(8,'--')]:axs[1,0].plot(S['t'][:,0],S['Pi_advection'][:,k==K].ravel(),ls,label=f'{label}, K={K}',alpha=.7)
 axs[1,1].plot(a['t'],a['lowpass_skewness'],label=label,alpha=.7)
for ax in axs.flat:ax.set_xlabel('Model time t');ax.legend(fontsize=7)
axs[0,0].set(title='Total kinetic energy',ylabel='E')
axs[0,1].set(title='Absolute box-scale energy, |k| ≤ 2',ylabel='E_low')
axs[1,0].set(title='Native advective flux: negative = upscale',ylabel='Π(K)');axs[1,0].axhline(0,color='k',lw=.5)
axs[1,1].set(title='Low-pass vorticity skewness, |k| ≤ 2',ylabel='Skewness');axs[1,1].axhline(0,color='k',lw=.5)
fig.suptitle('Ordinary viscosity only · no drag · phase-swept diagnostics')
fig.savefig(out/'diagnostics.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(2,1,figsize=(10,8),layout='constrained')
for case,label in zip(cases,labels):
 S,k,late=data[case][4],data[case][5],data[case][8]
 axs[0].loglog(k[1:],S['energy'][late,1:].mean(0),label=label)
 axs[1].plot(k,S['Pi_advection'][late].mean(0),label=label)
axs[0].set(ylabel='Shell kinetic energy',xlim=(1,128));axs[1].set(ylabel='Π(K)',xlim=(0,24));axs[1].axhline(0,color='k',lw=.5)
for ax in axs:ax.set_xlabel('Wavenumber');ax.axvspan(12,16,color='.8',alpha=.4);ax.legend()
fig.suptitle('Late window, t ≥ 100 · complete pulse-phase cycles')
fig.savefig(out/'spectra_flux.png',dpi=150);plt.close(fig)
def field(case,row):
 p,c=data[case][:2];q=np.fromfile(p/row['file'],dtype=np.float32).reshape(c['grid'],order='F');assert np.isfinite(q).all();return q
assert all(np.array_equal(field(c,data[c][6][0]),field('forced',data['forced'][6][0])) for c in cases)
limit=max(float(abs(field(c,row)).max()) for c in cases for row in data[c][6])
norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
fig,axs=plt.subplots(2,2,figsize=(12,11),layout='constrained');images=[]
for ax,c,label in zip(axs.flat,cases,labels):
 images.append(ax.imshow(field(c,data[c][6][0]).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),norm=norm,cmap='RdBu_r',interpolation='nearest'));ax.set(title=label,xlabel='x',ylabel='y')
fig.colorbar(images[0],ax=axs.ravel().tolist(),shrink=.8,label='Vertical vorticity ζ · fixed shared asinh scale')
text=fig.suptitle('');fig.supxlabel('Viscosity only · no drag · one horizontal 2D plane\nIdentical initial turbulence · positive ζ is counterclockwise · dimensionless model time',fontsize=10)
movie=out/'viscosity_only_four_controls.mp4'
writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-b:v','2200k','-maxrate','2400k','-bufsize','4800k','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=100):
 for i,t in enumerate(data['forced'][7]):
  for im,c in zip(images,cases):im.set_data(field(c,data[c][6][i]).T)
  text.set_text(f'Viscosity-only comparison · t = {t:.1f} / 200');writer.grab_frame()
  if i==400:fig.savefig(out/'poster.png',dpi=100)
plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==401
subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
report['movie']={'file':movie.name,'frames':count,'duration_seconds':duration,'full_decode_passed':True,'sampling':'fields every 0.5 model time, alternating post-pulse and mid-cycle; inherited config movie_sampling label is inaccurate'}
report['color_scale']={'type':'symmetric asinh','linear_width':.1,'limits':[-limit,limit],'fixed_shared':True,'clipping':False}
report['assets']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print('MOVIE_AND_DIAGNOSTICS_COMPLETE',flush=True)

"""Unity-only analysis/rendering of the completed matched continuation."""
from pathlib import Path
import csv,hashlib,json,subprocess,sys,tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
cases=['forced','unforced'];data={};summaries={}
for case in cases:
    p=src/case;c=tomllib.loads((p/'config.toml').read_text());r=tomllib.loads((p/'review.toml').read_text())
    assert r['status']=='COMPLETED', (case,r['status'])
    a=np.genfromtxt(p/'statistics.csv',delimiter=',',names=True)
    s=np.genfromtxt(p/'spectra.csv',delimiter=',',names=True)
    frames=list(csv.DictReader((p/'frames.csv').open()))
    ts=np.array([float(x['growth_times']) for x in frames])
    assert len(frames)==901 and np.isclose(ts[0],120) and np.isclose(ts[-1],1200)
    assert np.allclose(np.diff(ts),1.2)
    k=np.unique(s['k']);nt=len(s)//len(k)
    S={n:s[n].reshape(nt,len(k)) for n in s.dtype.names}
    assert np.all(S['k']==k[None,:])
    assert abs(S['T_injection'][:,k<6]).max()<1e-10
    windows={}
    for lo,hi in [(120,240),(240,600),(600,1200)]:
        sel=(S['growth_times'][:,0]>=lo-1e-7)&(S['growth_times'][:,0]<hi-1e-7)
        assert sel.any()
        windows[f'{lo}_{hi}']={
            'mean_flux_K1_to_K5':[float(S['Pi_advection'][sel][:,k==K].mean()) for K in range(1,6)],
            'fraction_negative_flux_K1_to_K5':[float((S['Pi_advection'][sel][:,k==K]<0).mean()) for K in range(1,6)],
            'mean_injection_work':float(S['T_injection'][sel].sum(1).mean()),
            'mean_drag_work':float(S['T_drag'][sel].sum(1).mean()),
            'mean_viscous_work':float(S['T_viscosity'][sel].sum(1).mean()),
            'mean_energy_spectrum':S['energy'][sel].mean(0).tolist(),
            'mean_flux':S['Pi_advection'][sel].mean(0).tolist()}
    echange=float(a['energy'][-1]-a['energy'][0]);budget=sum(r['energy_increments'])
    assert abs(echange-budget)<1e-8*float(a['energy'][0])
    summaries[case]={'config':c,'review':r,'windows':windows,'shell_upper_edges':k.tolist(),
        'energy_initial_final':[float(a['energy'][0]),float(a['energy'][-1])],
        'low_k_energy_initial_final':[float(a['energy'][i]*a['low_k_energy_fraction'][i]) for i in [0,-1]],
        'energy_storage_minus_integrated_terms':echange-budget}
    data[case]=(p,c,r,a,S,k,frames,ts)
assert summaries['forced']['config']['source_checkpoint_sha256']==summaries['unforced']['config']['source_checkpoint_sha256']

plt.rcParams.update({'font.size':10})
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for case in cases:
    p,c,r,a,S,k,frames,ts=data[case]
    axs[0,0].semilogy(a['growth_times'],a['energy'],label=case)
    axs[0,1].semilogy(a['growth_times'],a['energy']*a['low_k_energy_fraction'],label=case)
    axs[1,0].plot(a['growth_times'],a['integral_length'],label=case)
    for K,ls in [(2,'-'),(5,'--')]:
        # Fixed-width blocks improve visibility while retaining raw tabular data.
        t=S['growth_times'][:,0];f=S['Pi_advection'][:,k==K].ravel();step=20
        tb=[t[i:i+step].mean() for i in range(0,len(t),step)]
        fb=[f[i:i+step].mean() for i in range(0,len(t),step)]
        axs[1,1].plot(tb,fb,ls,label=f'{case}, K={K}')
axs[0,0].set(ylabel='Kinetic energy',title='Total energy')
axs[0,1].set(ylabel='Kinetic energy at |k| ≤ 2',title='Absolute large-scale energy, not just its fraction')
axs[1,0].set(ylabel='Integral wavelength',title='2π × energy-weighted mean 1/|k|')
axs[1,1].set(ylabel='Advection energy flux Π(K)',title='12-growth-time block means; negative = upscale')
axs[1,1].axhline(0,c='k',lw=.8)
for ax in axs.flat:ax.set_xlabel('Reference growth times σt');ax.legend()
fig.suptitle('Matched 2D forcing-on/off continuations from the same turbulent state')
fig.savefig(out/'energy_transfer_history.png',dpi=150);plt.close(fig)

fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for col,case in enumerate(cases):
    k=data[case][5]
    for name,w in summaries[case]['windows'].items():
        axs[0,col].loglog(k[1:],np.array(w['mean_energy_spectrum'])[1:],label=name.replace('_','–'))
        axs[1,col].plot(k,w['mean_flux'],label=name.replace('_','–'))
    axs[0,col].set(title=case,ylabel='Shell kinetic energy',xlim=(1,64))
    axs[1,col].set(title='Negative flux = upscale',ylabel='Π(K)',xlim=(0,12))
    axs[1,col].axhline(0,c='k',lw=.8)
    for ax in axs[:,col]:ax.axvspan(6,8,color='.8',alpha=.5);ax.set_xlabel('Wavenumber / shell upper edge');ax.legend()
fig.suptitle('Spectra and nonlinear energy flux · shaded region = forcing band')
fig.savefig(out/'spectra_flux.png',dpi=150);plt.close(fig)

def field(case,row):
    p,c,*_=data[case];q=np.fromfile(p/row['file'],dtype=np.float32).reshape(c['grid'],order='F')
    assert np.isfinite(q).all();return q
assert np.array_equal(field('forced',data['forced'][6][0]),field('unforced',data['unforced'][6][0]))
limit=max(float(abs(field(case,row)).max()) for case in cases for row in data[case][6])
norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ffmpeg
fig,axs=plt.subplots(1,2,figsize=(12,6.5),layout='constrained');images=[]
for ax,case in zip(axs,cases):
    im=ax.imshow(field(case,data[case][6][0]).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),norm=norm,cmap='RdBu_r',interpolation='nearest')
    ax.set(title='Forcing ON' if case=='forced' else 'Forcing OFF',xlabel='x',ylabel='y');images.append(im)
fig.colorbar(images[0],ax=axs,shrink=.8,label='Vertical vorticity ζ · fixed shared asinh scale')
title=fig.suptitle('');fig.supxlabel('Single horizontal 2D plane · identical turbulent start at σt=120\nPositive ζ is counterclockwise; σ=0.6 is a common reference clock, not Earth/Jupiter days',fontsize=10)
movie=out/'forced_unforced_120_1200.mp4';writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-crf','25','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=120):
    for i,t in enumerate(data['forced'][7]):
        assert np.isclose(t,data['unforced'][7][i])
        for case,im in zip(cases,images):im.set_data(field(case,data[case][6][i]).T)
        title.set_text(f'2D matched continuations · σt = {t:.1f} / 1200')
        writer.grab_frame()
        if i==len(data['forced'][7])-1:fig.savefig(out/'poster.png',dpi=120)
plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==901
subprocess.run([ffmpeg,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
report={'status':'COMPLETE_REQUIRES_SCIENTIFIC_REVIEW','compute_location':'Unity','cases':summaries,
    'movie':{'file':movie.name,'frames':count,'fps':12,'duration_seconds':duration,'full_decode_passed':True},
    'color_scale':{'type':'symmetric asinh','linear_width':.1,'limits':[-limit,limit],'fixed_shared':True,'clipping':False},
    'limitations':['This is a mature-state forcing switch at sigma*t=120, not forcing-off from the original weak noise.',
        'Only a short range lies below the forcing band; negative flux does not establish a broad inertial range.',
        'A growing low-k energy fraction can reflect selective decay; absolute energy and native advective flux are reported.',
        'Native numerical-budget closure is not resolution convergence or a sign-selection mechanism.',
        'The unforced reference clock is defined using the forced-case sigma; no physical-day mapping.'],
    'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}}
(out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':'COMPLETE','frames':count,'duration_seconds':duration,'movie_bytes':movie.stat().st_size}),flush=True)

"""Completed stress-suite analysis and rendering; run only on Unity."""
from pathlib import Path
import csv,json,sys,tomllib,subprocess,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

root,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
cases=['decay','symmetric','biased_plus','biased_minus']
labels={'decay':'Ordinary decay','symmetric':'Symmetric stress · f = +100','biased_plus':'Sign-biased stress · f = +100','biased_minus':'Sign-biased stress · f = −100 · reflected IC'}
frames={};scalars={};configs={}
for case in cases:
    p=root/case;c=tomllib.loads((p/'config.toml').read_text());s=tomllib.loads((p/'status.toml').read_text())
    assert s['status']=='COMPLETED' and abs(s['t']-400)<1e-7
    assert c['case']==case and c['grid']==[256,256] and c['viscosity']==.0001 and c['cutoff']==16 and c['beta']==0
    assert c['nu0']==(0 if case=='decay' else -.0002)
    assert c['nu1']==(-.001 if case.startswith('biased') else 0)
    assert c['nu2']==(0 if case=='decay' else .01)
    rows=list(csv.DictReader((p/'frames.csv').open()));assert [float(x['t']) for x in rows]==list(range(401))
    frames[case]=[p/x['file'] for x in rows]
    a=np.genfromtxt(p/'scalars.csv',delimiter=',',names=True);assert np.array_equal(a['t'],np.arange(401))
    for name in a.dtype.names:assert np.isfinite(a[name]).all(),(case,name)
    for part in ['linear','quadratic','cubic']:
        assert np.allclose(a[part+'_energy_rate'],a[part+'_energy_identity'],rtol=1e-9,atol=1e-12)
    assert (a['cubic_energy_rate']<=1e-12).all()
    scalars[case]=a;configs[case]=c
assert len({c['initial_sha256'] for c in configs.values()})==1
def field(case,i):
    q=np.fromfile(frames[case][i],dtype=np.float32);assert q.size==256**2 and np.isfinite(q).all()
    return q.reshape((256,256),order='F')
assert np.array_equal(field('decay',0),field('symmetric',0))
assert np.array_equal(field('decay',0),field('biased_plus',0))
idx=np.r_[0,np.arange(255,0,-1)]
assert np.array_equal(field('biased_minus',0),-field('biased_plus',0)[idx,:])
limit=max(float(np.abs(field(c,i)).max()) for c in cases for i in range(401))
norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
movies={}
for selected,stem in [(cases,'stress_controls_256_0_400'),(['biased_plus'],'stress_plus_256_0_400')]:
    rows,cols=(2,2) if len(selected)==4 else (1,1)
    fig,axs=plt.subplots(rows,cols,figsize=(12,11) if rows==2 else (7.6,6.8),squeeze=False,layout='constrained')
    images=[]
    for ax,case in zip(axs.flat,selected):
        images.append(ax.imshow(field(case,0).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),cmap='RdBu_r',norm=norm,interpolation='nearest'))
        ax.set(title=labels[case],xlabel='x',ylabel='y')
    fig.colorbar(images[0],ax=list(axs.flat),shrink=.8,label='Vertical vorticity ζ · fixed shared asinh scale')
    title=fig.suptitle('');fig.supxlabel('256² · horizontal 2D model · dimensionless time\nCyclonic: positive ζ for +f, negative ζ for −f. No tracer or depth layers.',fontsize=10)
    movie=out/(stem+'.mp4');writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-b:v','1800k','-maxrate','2200k','-bufsize','4400k','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
    with writer.saving(fig,str(movie),dpi=100):
        for i in range(401):
            for im,case in zip(images,selected):im.set_data(field(case,i).T)
            title.set_text(f'Projected nonlinear-stress hypothesis · t = {i} / 400');writer.grab_frame()
            if i==400:fig.savefig(out/(stem+'.png'),dpi=100)
    plt.close(fig)
    count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==401 and abs(duration-401/12)<.1
    subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
    assert movie.stat().st_size<25*1024**2
    movies[stem]={'frames':count,'duration_seconds':duration,'full_decode_passed':True}
    print('MEDIA_VERIFIED',stem,flush=True)
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for ax,key,title in zip(axs.flat,['energy_k_le_2','skewness','filtered_rossby_rms','filtered_rossby_max'],['Absolute energy in |k| ≤ 2','Vorticity skewness','Filtered RMS |Z/f|','Filtered maximum |Z/f|']):
    for case in cases:ax.plot(scalars[case]['t'],scalars[case][key],label=labels[case])
    ax.set(xlabel='Dimensionless time',title=title);ax.grid(alpha=.2)
axs[0,0].legend(fontsize=8);axs[1,1].axhline(.2,color='gray',ls='--',label='0.2 reference, not a validity theorem')
fig.suptitle('Stress-model experiment · descriptive diagnostics, not proof of an inverse cascade')
fig.savefig(out/'diagnostics.png',dpi=150);plt.close(fig)
summary={}
for case,a in scalars.items():
    summary[case]={'final':{name:float(a[name][-1]) for name in a.dtype.names},
        'maximum_filtered_abs_rossby':float(a['filtered_rossby_max'].max()),
        'maximum_fraction_abs_rossby_gt_02':float(a['filtered_fraction_abs_ro_gt_02'].max()),
        'maximum_filtered_strain_over_f':float(a['filtered_strain_max_over_f'].max()),
        'maximum_highk_enstrophy_fraction':float(a['highk_enstrophy_fraction'].max())}
    (out/(case+'_scalars.csv')).write_bytes((root/case/'scalars.csv').read_bytes())
public={'status':'COMPLETED_MEDIA_VERIFIED','grid':[256,256],'domain':[float(2*np.pi)]*2,'time_range':[0,400],
 'time_units':'dimensionless','compute_location':'Unity','viscosity':.0001,'cutoff':16,'beta':0,
 'cases':{k:{n:c[n] for n in ['nu0','nu1','nu2','f','reflected_initial_condition']} for k,c in configs.items()},
 'movies':movies,'color_scale':{'type':'symmetric asinh','width':.1,'limits':[-limit,limit],'fixed_across_cases_and_time':True,'clipping':False},
 'summary':summary,'source_work_identity_passed':True,'matched_or_reflected_initial_fields':True,
 'limitations':['Constitutive hypothesis, not resolved convection or a closure derived from Breeze.',
 'Stress acts directly on low modes; large-scale growth and morphology do not alone demonstrate an inverse cascade.',
 'Only one resolution and one initial realization; reflected minus-f case is a symmetry check, not an independent sample.',
 'Small Rossby number is only one requirement of the formal expansion; scale separation and rapid adjustment are assumptions.',
 'Source-work checks are algebraic identities, not an integrated budget of all numerical tendencies.']}
(out/'report.json').write_text(json.dumps(public,indent=2,allow_nan=False)+'\n')
(out/'transport.json').write_text(json.dumps({'source_root':str(root),'configs':configs,'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}},indent=2)+'\n')
print('PUBLICATION_READY',flush=True)

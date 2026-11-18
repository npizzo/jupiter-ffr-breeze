"""Unity-only matched formation movie and final diagnostics, after all segments succeed."""
from pathlib import Path
import sys,csv,json,tomllib,subprocess,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'feedback_maintenance_v1'))
import analyze_existing as old
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

def main(root,out):
    out.mkdir(parents=True,exist_ok=False)
    original,_,cfg,_=old.inventory(Path('data/viscous2d-feedback512-v1-1925368'),'biased_plus')
    frames={};scalars={};configs=[]
    for end in (80,160,240,320,400):
        d=root/f'to_{end:03d}';c=tomllib.loads((d/'config.toml').read_text());status=tomllib.loads((d/'status.toml').read_text())
        assert status['status']=='COMPLETED' and status['t']==end
        assert c['initial_sha256']==cfg['initial_sha256'] and c['momentum_advection'] is False
        assert c['a']==.0002 and c['b']==.002 and c['c']==.0002 and c['viscosity']==.0001 and c['beta']==0
        for r in csv.DictReader((d/'frames.csv').open()):
            t=float(r['t']);p=d/r['file']
            if t in frames:assert p.read_bytes()==frames[t].read_bytes()
            frames[t]=p
        for r in csv.DictReader((d/'scalars.csv').open()):
            row={k:float(v) for k,v in r.items()};scalars[row['t']]=row
        configs.append(c)
    times=np.arange(0,400.5,.5);assert sorted(frames)==list(times)
    assert np.array_equal(old.load_field(original[0.]),old.load_field(frames[0.]))
    st=np.array(sorted(scalars));assert np.array_equal(st,np.arange(0,400.25,.25))
    rows=[scalars[t] for t in st]
    assert all(np.isfinite(list(r.values())).all() for r in rows)
    assert max(abs(r['advection_lowk_energy_rate']) for r in rows)<1e-14
    potential=np.array([r['potential'] for r in rows]);increase=float(np.diff(potential).max())
    assert increase<1e-9*max(float(np.abs(potential).max()),1.),increase
    integral=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
    work={k:float(integral([r[k+'_lowk_energy_rate'] for r in rows],st)) for k in ('advection','linear','quadratic','cubic','viscous')}
    change=rows[-1]['energy_k_le_2']-rows[0]['energy_k_le_2']
    limit=max(float(np.abs(old.load_field(p)).max()) for files in (original,frames) for p in files.values())
    norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
    ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
    fig,axs=plt.subplots(1,2,figsize=(13.6,6.8),layout='constrained');images=[]
    for ax,label,files in zip(axs,('Original: advection on','Matched control: advection off'),(original,frames)):
        images.append(ax.imshow(old.load_field(files[0.]).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),cmap='RdBu_r',norm=norm,interpolation='nearest'))
        ax.set(title=label,xlabel='x',ylabel='y')
    fig.colorbar(images[0],ax=list(axs),label='Vorticity ζ · fixed shared asinh scale',shrink=.8)
    title=fig.suptitle('');fig.supxlabel('512² · same original initial condition and polynomial coefficients · f=+1 · beta=0\nHorizontal 2D: one plane, no tracer · positive vorticity is cyclonic',fontsize=10)
    movie=out/'formation_advection_on_off_0_400.mp4'
    writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-b:v','1600k','-maxrate','1800k','-bufsize','3600k','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
    with writer.saving(fig,str(movie),dpi=100):
        for i,t in enumerate(times):
            images[0].set_data(old.load_field(original[t]).T);images[1].set_data(old.load_field(frames[t]).T)
            title.set_text(f'Polynomial-model formation · t = {t:g}')
            writer.grab_frame()
            if i==len(times)-1:fig.savefig(out/'formation_advection_on_off_0_400.png',dpi=100)
    plt.close(fig)
    count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==801 and abs(duration-801/12)<.1
    subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
    assert movie.stat().st_size<25*1024**2
    report={'status':'COMPLETED_MEDIA_VERIFIED','compute_location':'Unity','time_range':[0,400],'intervention':'Only native momentum advection disabled',
       'no_advection':{'initial':rows[0],'final':rows[-1],'integrated_lowmode_work':work,'energy_change':change,'sampling_residual':change-sum(work.values()),'maximum_potential_increment':increase},
       'movie':{'file':movie.name,'frames':count,'seconds':duration,'decode_pass':True,'layout':'1x2 synchronized horizontal vorticity','color_scale':{'type':'asinh','width':.1,'limits':[-limit,limit],'fixed':True}},
       'limitations':['One initial realization; formation control, not a validated convection closure.','Cumulative work is sampled at 0.25, not an exact every-stage budget.']}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
    (out/'transport.json').write_text(json.dumps({'assets':assets,'configs':configs},indent=2));print('FORMATION_MOVIE_READY',flush=True)

if __name__=='__main__':main(*map(Path,sys.argv[1:]))

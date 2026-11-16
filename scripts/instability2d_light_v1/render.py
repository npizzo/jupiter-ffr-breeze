"""Render saved 2D vorticity on Unity; no inferred vertical levels or tracer."""
from pathlib import Path
import csv, hashlib, json, subprocess, sys, tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
config=tomllib.loads((src/'config.toml').read_text())
review=tomllib.loads((src/'review.toml').read_text())
rows=list(csv.DictReader((src/'frames.csv').open()))
nx,ny=config['grid'];times=np.array([float(r['growth_times']) for r in rows])
assert len(rows)==101 and np.all(np.diff(times)>0)
assert np.isclose(times[0],0) and np.isclose(times[-1],120)
def field(row):
    q=np.fromfile(src/row['file'],dtype=np.float32).reshape((nx,ny),order='F')
    assert np.all(np.isfinite(q))
    return q
limit=max(float(np.max(np.abs(field(r)))) for r in rows)
assert limit>0
ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
plt.rcParams['animation.ffmpeg_path']=ffmpeg
plt.rcParams.update({'font.size':11})
fig,ax=plt.subplots(figsize=(8,7.6),layout='constrained')
norm=AsinhNorm(linear_width=0.1,vmin=-limit,vmax=limit)
im=ax.imshow(field(rows[0]).T,origin='lower',extent=[0,2*np.pi,0,2*np.pi],
             cmap='RdBu_r',norm=norm,interpolation='nearest',aspect='equal')
ax.set(xlabel='x (dimensionless)',ylabel='y (dimensionless)')
cb=fig.colorbar(im,ax=ax,shrink=.86,label=r'Vertical vorticity $\zeta$ (dimensionless)')
ticks=[-limit]+[x for x in [-10,-3,-1,-.3,0,.3,1,3,10] if abs(x)<limit]+[limit]
cb.set_ticks(ticks);cb.set_ticklabels([f'{x:.2g}' for x in ticks])
fig.suptitle('2D instability-forced Oceananigans pilot\n128 × 128 · single horizontal plane · no beta',fontsize=14)
title=ax.set_title('')
fig.supxlabel('Fixed symmetric asinh scale · positive vorticity is counterclockwise\nTime is measured in forcing growth times, not Earth/Jupiter days',fontsize=10)
movie=out/'vorticity_0_120_growth_times.mp4'
writer=FFMpegWriter(fps=8,codec='libx264',extra_args=['-crf','26','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=120):
    for i,row in enumerate(rows):
        im.set_data(field(row).T)
        title.set_text(f"Growth times: {float(row['growth_times']):6.1f} / 120     t = {float(row['t']):6.1f}")
        writer.grab_frame()
        if i==len(rows)-1:fig.savefig(out/'poster.png',dpi=120)
plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie))
assert count==len(rows)
subprocess.run([ffmpeg,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
metadata={'status':'COMPLETE','compute_location':'Unity','simulation_job':'1810727',
    'config':config,'review':review,'movie':{'file':movie.name,'frames':count,'duration_seconds':duration,
    'fps':8,'growth_times':[float(times[0]),float(times[-1])],'full_decode_passed':True},
    'color_scale':{'type':'symmetric asinh','linear_width':0.1,'limits':[-limit,limit],
                   'fixed_across_all_frames':True,'clipping':False,'colormap':'RdBu_r'},
    'scope':'Single-plane 2D pilot; no depth or tracer dimension exists. Not a three-depth convection run.',
    'limitation':'Native budgets close, but numerical-advection metric 10.53% exceeds 10% attribution cutoff; no physical sign-selection mechanism established.',
    'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
              for p in out.iterdir() if p.is_file()}}
(out/'metadata.json').write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
print(json.dumps({'status':'PASS','frames':count,'duration_seconds':duration,
                  'movie_bytes':movie.stat().st_size,'limits':[-limit,limit]}),flush=True)

"""Render the frozen t=0–160 half-wavelength prefix on Unity."""
from pathlib import Path
import csv,hashlib,json,subprocess,sys,tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg
original,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
frames={};segments=[];duplicates=[]
for p,lo,hi,dt in [(original/'to_080',0,80,.00125),(original/'to_160',80,160,.00125)]:
    c=tomllib.loads((p/'config.toml').read_text())
    status=tomllib.loads((p/'status.toml').read_text());assert status['status']=='COMPLETED' and abs(status['t']-hi)<1e-7
    for k,v in {'case':'halfscale','grid':[512,512],'domain':[2*np.pi,2*np.pi],'f':1.,'beta':0.,'dt':dt,'viscosity':.0001,'drag':0.,'instability_coefficient':.000175,'forcing_band':[24.,32.],'allow_zonal_forcing':False,'frame_interval':.5}.items():
        assert c[k]==v,(k,c[k],v)
    selected=[]
    for row in csv.DictReader((p/'frames.csv').open()):
        it,t=int(row['iteration']),float(row['t']);idx=round(t*2)
        if idx<2*lo or idx>2*hi:continue
        assert abs(t-it*dt)<1e-7 and abs(t-idx*.5)<1e-7
        path=p/row['file'];assert path.stat().st_size==512*512*4
        if idx in frames:
            assert path.read_bytes()==frames[idx][1].read_bytes(),'Restart join differs'
            duplicates.append(idx*.5)
        frames[idx]=(t,path);selected.append(idx)
    assert selected==list(range(lo*2,hi*2+1)),(lo,hi,len(selected))
    segments.append({'directory':str(p),'published_range':[lo,hi],'config':c})
assert sorted(frames)==list(range(321)) and duplicates==[80.]
assert len({s['config']['initial_sha256'] for s in segments})==1
paths=[frames[i][1] for i in range(321)]
def field(i):
    q=np.fromfile(paths[i],dtype=np.float32);assert q.size==512*512 and np.isfinite(q).all()
    return q.reshape((512,512),order='F')
limit=max(float(np.abs(field(i)).max()) for i in range(321));assert limit>0
norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
fig,ax=plt.subplots(figsize=(8.4,7.4),layout='constrained')
im=ax.imshow(field(0).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),cmap='RdBu_r',norm=norm,interpolation='nearest')
ax.set(xlabel='x',ylabel='y');title=ax.set_title('')
fig.colorbar(im,ax=ax,shrink=.83,label='Vertical vorticity ζ · fixed asinh scale')
fig.supxlabel('512 × 512 · a = 0.000175 · band 24–32 · ordinary viscosity only · no drag\nOne horizontal 2D plane · positive ζ is counterclockwise · dimensionless time',fontsize=9)
stem='halfscale_512_0_160';movie=out/(stem+'.mp4')
writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-b:v','2200k','-maxrate','2400k','-bufsize','4800k','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
with writer.saving(fig,str(movie),dpi=100):
    for i in range(321):
        im.set_data(field(i).T);title.set_text(f'Half-wavelength instability forcing · t = {i*.5:.1f} / 160')
        writer.grab_frame()
        if i==320:fig.savefig(out/(stem+'.png'),dpi=100)
        if i%100==0:print('FRAME',i,'of 321',flush=True)
plt.close(fig)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==321 and abs(duration-321/12)<.1
subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
assert movie.stat().st_size<25*1024**2
public={'status':'MEDIA_VERIFIED','published_time_range':[0,160],'time_units':'dimensionless model time','frames':count,'fps':12,'duration_seconds':duration,'full_decode_passed':True,'duplicate_restart_times_verified_equal':duplicates,'join_by':'physical time; duplicate t=80 restart frame shown once','color_scale':{'type':'symmetric asinh','linear_width':.1,'limits':[-limit,limit],'fixed_across_time':True,'clipping':False},'parameters':{k:segments[0]['config'][k] for k in ['grid','domain','f','beta','viscosity','drag','instability_coefficient','forcing_band','allow_zonal_forcing','frame_interval']},'timestep_segments':[{'time_range':[0,160],'dt':.00125}],'notes':['Single horizontal 2D plane; no depth layers or passive tracer.','All 321 saved times, without interpolation; restart joins shown once.','Forcing wavelength is halved; coefficient adjusted to match continuum net growth after viscosity to the doubled-forcing reference. Initial velocities are unchanged.','This is a frozen publication through t=160, not a claim of completion through t=400.','No inverse-energy-transfer or vortex-sign mechanism inferred from morphology alone.']}
(out/'report.json').write_text(json.dumps(public,indent=2)+'\n')
private={'segments':segments,'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.suffix in ['.mp4','.png','.json']}}
(out/'transport.json').write_text(json.dumps(private,indent=2)+'\n')
print('MEDIA_COMPLETE',count,duration,flush=True)

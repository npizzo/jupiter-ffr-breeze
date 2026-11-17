"""Render the completed controlled polynomial-feedback suite on Unity."""
from pathlib import Path
import csv,json,subprocess,sys,tomllib,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg
src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
cases=['linear_only','cubic_control','biased_plus','biased_minus']
labels={'linear_only':'Existing linear instability only · f = +1','cubic_control':'Linear instability + cubic saturation · f = +1','biased_plus':'Linear + quadratic + cubic · f = +1','biased_minus':'Reflected initial flow · same closure · f = −1'}
data={};private={}
for case in cases:
    if case=='linear_only':
        pairs=[(n,Path('data/viscous2d-long512-v1-1848831/forced')/f'to_{n:03d}') for n in [40,80,120,160,200]]+[(400,Path('data/viscous2d-fast512-v1-1872263/forced'))]
    else:pairs=[(n,src/case/f'to_{n:03d}') for n in [80,160,240,320,400]]
    frames={};segments=[]
    for end,p in pairs:
        c=tomllib.loads((p/'config.toml').read_text());s=tomllib.loads((p/'status.toml').read_text())
        assert s['status']=='COMPLETED' and abs(s['t']-end)<1e-7
        assert c['grid']==[512,512] and c['viscosity']==.0001 and c['beta']==0
        if case!='linear_only':
            assert c['case']==case and c['a']==.0002 and c['c']==.0002
            assert c['b']==(0 if case=='cubic_control' else .002)
            assert c['f']==(-1 if case=='biased_minus' else 1)
        else:assert c['instability_coefficient']==.0002 and c['drag']==0
        for row in csv.DictReader((p/'frames.csv').open()):
            t=float(row['t']);ix=round(t*2);assert abs(ix*.5-t)<1e-7
            path=p/row['file']
            if ix in frames:assert path.read_bytes()==frames[ix].read_bytes()
            frames[ix]=path
        segments.append({'directory':str(p),'config':c})
    assert sorted(frames)==list(range(801));data[case]=[frames[i] for i in range(801)];private[case]=segments
assert len({s['config']['initial_sha256'] for seg in private.values() for s in seg})==1
def field(case,i):
    q=np.fromfile(data[case][i],dtype=np.float32);assert q.size==512*512 and np.isfinite(q).all()
    return q.reshape((512,512),order='F')
assert np.array_equal(field('linear_only',0),field('cubic_control',0))
assert np.array_equal(field('linear_only',0),field('biased_plus',0))
# C-grid corner curl under x reflection: minus original at reflected x-face indices.
faceidx=np.r_[0,np.arange(511,0,-1)]
assert np.array_equal(field('biased_minus',0),-field('biased_plus',0)[faceidx,:])
limit=max(float(np.abs(field(c,i)).max()) for c in cases for i in range(801))
norm=AsinhNorm(linear_width=.1,vmin=-limit,vmax=limit)
ff=imageio_ffmpeg.get_ffmpeg_exe();plt.rcParams['animation.ffmpeg_path']=ff
movies={}
for selected,stem in [(cases,'feedback_controls_512_0_400'),(['biased_plus'],'feedback_plus_512_0_400')]:
    rows,cols=(2,2) if len(selected)==4 else (1,1)
    fig,axs=plt.subplots(rows,cols,figsize=(13,12) if rows==2 else (7.6,6.8),squeeze=False,layout='constrained')
    images=[]
    for ax,case in zip(axs.flat,selected):
        images.append(ax.imshow(field(case,0).T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),cmap='RdBu_r',norm=norm,interpolation='nearest'))
        ax.set(title=labels[case],xlabel='x',ylabel='y')
    fig.colorbar(images[0],ax=list(axs.flat),shrink=.8,label='Vertical vorticity ζ · fixed shared asinh scale')
    title=fig.suptitle('');fig.supxlabel('512² · one horizontal 2D plane · dimensionless time\nSigned vorticity shown: cyclonic is positive for f=+1 and negative for f=−1',fontsize=10)
    movie=out/(stem+'.mp4');writer=FFMpegWriter(fps=12,codec='libx264',extra_args=['-b:v','2200k','-maxrate','2400k','-bufsize','4800k','-threads','1','-pix_fmt','yuv420p','-movflags','+faststart'])
    with writer.saving(fig,str(movie),dpi=100):
        for i in range(801):
            for im,case in zip(images,selected):im.set_data(field(case,i).T)
            title.set_text(f'Polynomial-feedback hypothesis · t = {i*.5:.1f} / 400');writer.grab_frame()
            if i==800:fig.savefig(out/(stem+'.png'),dpi=100)
    plt.close(fig)
    count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==801 and abs(duration-66.75)<.1
    subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
    assert movie.stat().st_size<25*1024**2
    movies[stem]={'frames':count,'duration_seconds':duration,'full_decode_passed':True}
    print('VERIFIED',movie.name,flush=True)
public={'status':'MEDIA_VERIFIED','time_range':[0,400],'time_units':'dimensionless model time','compute_location':'Unity','frame_interval':.5,'movies':movies,'initial_fields_matched_or_reflected':True,'color_scale':{'type':'symmetric asinh','width':.1,'limits':[-limit,limit],'fixed_across_cases_and_time':True,'clipping':False},'parameters':{'grid':[512,512],'a':.0002,'b':.002,'c':.0002,'viscosity':.0001,'forcing_band':[12,16],'beta':0,'dt_ceiling':.0025},'notes':['Cubic control has b=0. Existing linear-only control has b=c=0.','Minus-f case starts from x-reflected initial velocity; it is a symmetry test, not an independent realization.','Source is mean-removed b*f*zeta^2-c*zeta^3. It is imposed phenomenology, not a derived convection or resonance model.','Scalar records distinguish large-scale energy and direct polynomial energy input; morphology alone does not establish inverse cascade.','No tracer or vertical layers exist in these horizontal 2D cases.']}
(out/'report.json').write_text(json.dumps(public,indent=2)+'\n')
(out/'transport.json').write_text(json.dumps({'segments':private,'assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.suffix in ['.mp4','.png','.json']}},indent=2)+'\n')
print('MEDIA_COMPLETE',flush=True)

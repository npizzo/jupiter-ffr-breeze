"""Render three real vorticity planes, with fixed scales across time and matched cases."""
import argparse
import json
from pathlib import Path
import tomllib
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter


def read(path, f, level=None):
    with h5py.File(path) as d:
        t=float(d['t'][()]); heights=d['video_heights'][:]
        field=d['vorticity_slices'][:] if level is None else d['vorticity_slices'][level]
    if not np.isfinite(field).all() or not np.isfinite(t):
        raise ValueError(f'Nonfinite output: {path}')
    return t, heights, field / f


def collect(folder, preview):
    segments=sorted(folder.glob('segment_*'), key=lambda p: int(p.name.removeprefix('segment_'))) if not (folder/'metadata.toml').exists() else [folder]
    segments=[s for s in segments if (s/'metadata.toml').exists()]
    if not segments:
        raise ValueError(f'No metadata found in {folder}')
    files=[]; records=[]
    for segment in segments:
        m=tomllib.loads((segment/'metadata.toml').read_text()); records.append(m)
        if not preview and m['status'] != 'COMPLETE':
            raise ValueError('Incomplete run: use --preview to label a partial time series')
        files.extend(sorted(segment.glob('frame_*.jld2')))
    m=records[-1]; times=[]; limit=np.zeros(3)
    if not np.isfinite(m['f']) or m['f']==0:
        raise ValueError('Nonzero finite f required')
    for file in files:
        t,heights,zeta=read(file,m['f'])
        if zeta.shape != (3,m['ny'],m['nx']) or not np.array_equal(heights,m['video_heights_actual_m']):
            raise ValueError(f'Invalid three-depth fields: {file}')
        times.append(t)
        limit=np.maximum(limit,np.max(np.abs(zeta),axis=(1,2)))
    times=np.asarray(times)
    if len(times)<2 or not np.all(np.diff(times)>0):
        raise ValueError('Need at least two distinct, chronological frames')
    expected=np.arange(records[0].get("segment_start_seconds",0),times[-1]+m['output_seconds']/2,m['output_seconds'])
    if times.shape != expected.shape or not np.allclose(times,expected,rtol=0,atol=1e-6):
        raise ValueError('Missing frames: movies must include every saved frame from the initial state')
    if not preview and not np.isclose(times[-1],m['stop_seconds']):
        raise ValueError('Final output missing')
    return m,files,times,limit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('cases',nargs='+',type=Path)
    p.add_argument('--output',type=Path,required=True,help='Directory for three MP4 videos, posters and provenance')
    p.add_argument('--labels',nargs='+',required=True)
    p.add_argument('--preview',action='store_true')
    p.add_argument('--fps',type=int,default=12)
    p.add_argument('--crf',type=int,default=30)
    a=p.parse_args()
    if len(a.cases) not in (1,2) or len(a.labels)!=len(a.cases) or a.fps<1:
        p.error('Provide one or two cases, one label per case, and a positive frame rate')
    cases=[collect(case,a.preview) for case in a.cases]
    ref=cases[0]
    for m,files,times,limits in cases[1:]:
        for key in ('nx','ny','nz','lx','ly','height','video_heights_actual_m','f','heat_flux','cooling_time','cooling_model'):
            if m[key]!=ref[0][key]: raise ValueError(f'Matched-case mismatch: {key}')
        if not np.array_equal(times,ref[2]): raise ValueError('Matched cases require identical saved times')
    limits=np.maximum(np.max([c[3] for c in cases],axis=0),1e-12)
    a.output.mkdir(parents=True,exist_ok=True)
    manifest={'labels':a.labels,'cases':[c[0]['name'] for c in cases], 'preview':a.preview,
              'times_seconds':ref[2].tolist(),'fps':a.fps,'videos':[],
              'note':'Real saved vertical vorticity, normalized by central f; no temporal interpolation. Each depth uses one symmetric maximum over all times and cases, with no clipping. Heights are above model bottom.'}
    for level,height in enumerate(ref[0]['video_heights_actual_m']):
        target=ref[0]['video_heights_requested_m'][level]
        filename=f'vorticity_{target/1000:g}km.mp4'
        fig,axes=plt.subplots(1,len(cases),figsize=(6*len(cases),5),squeeze=False,layout='constrained')
        images=[]
        for ax,label,(m,files,_,_) in zip(axes.flat,a.labels,cases):
            _,_,field=read(files[0],m['f'],level)
            image=ax.imshow(field,origin='lower',interpolation='nearest',cmap='RdBu_r',
                vmin=-limits[level],vmax=limits[level],extent=(0,m['lx']/1000,-m['ly']/2000,m['ly']/2000))
            ax.set(xlabel='x [km]',ylabel='y [km]',title=f'{label}\nz = {height/1000:.2f} km above bottom')
            images.append(image)
        fig.colorbar(images[0],ax=list(axes.flat),shrink=.8,label='Vertical vorticity ζ/f₀ (fixed scale)')
        title=fig.suptitle('')
        writer=FFMpegWriter(fps=a.fps,codec='libx264',extra_args=['-crf',str(a.crf),'-pix_fmt','yuv420p','-movflags','+faststart'])
        with writer.saving(fig,str(a.output/filename),dpi=100):
            for i,t in enumerate(ref[2]):
                for im,(m,files,_,_) in zip(images,cases): im.set_data(read(files[i],m['f'],level)[2])
                title.set_text(f'{"Partial run | " if a.preview else ""}Vorticity evolution | Earth day {t/86400:.3f}')
                writer.grab_frame()
        poster=filename.replace('.mp4','.png');fig.savefig(a.output/poster,dpi=100);plt.close(fig)
        manifest['videos'].append({'file':filename,'poster':poster,'height_m':float(height),
                                   'requested_height_m':target,'color_limits':[-float(limits[level]),float(limits[level])]})
    (a.output/'videos.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Rendered three depths, {len(ref[2])} real frames each, to {a.output}')

if __name__=='__main__':main()

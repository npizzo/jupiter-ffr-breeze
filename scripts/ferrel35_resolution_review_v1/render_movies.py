"""Two matched cases on one timeline, three depth columns, with shared fixed scales."""
import argparse
import json
import hashlib
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import AsinhNorm, FuncNorm
from vorticity_depth_videos_ordered import collect

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('cases', nargs='+', type=Path)
p.add_argument('--labels', nargs='+', required=True)
p.add_argument('--title', required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--tracer', action='store_true')
p.add_argument('--fps', type=int, default=12)
p.add_argument('--crf', type=int, default=29)
p.add_argument('--vorticity-width', type=float, default=0.1)
p.add_argument('--case-index',type=int,choices=(0,1),default=None,help='Render one case using limits calculated from both cases')
a = p.parse_args()
if len(a.cases) not in (1,2) or len(a.labels)!=len(a.cases):
    p.error('Provide one or two cases and one label per case')
import imageio_ffmpeg
plt.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()
if a.fps < 1 or a.vorticity_width <= 0:
    p.error('Positive fps and vorticity width required')
cases = [collect(path, False) for path in a.cases]
m, _, times, _ = cases[0]
for other, _, other_times, _ in cases[1:]:
    for key in ('nz', 'lx', 'ly', 'height', 'cooling_model',
                'cooling_time', 'video_heights_actual_m', 'turbulence_seed',
                'turbulence_wavelength'):
        if m[key] != other[key]:
            raise ValueError(f'Incompatible matched comparison: {key}')
    if not np.array_equal(times, other_times):
        raise ValueError('Cases require identical saved times')
if len(cases)==2:
    assert cases[1][0]['nx']==2*cases[0][0]['nx']
    assert cases[1][0]['ny']==2*cases[0][0]['ny']
limits = np.maximum(np.max([case[3] for case in cases], axis=0), 1e-12)
crange = [0., 1.]
if a.tracer:
    for meta, files, _, _ in cases:
        for file in files:
            with h5py.File(file) as d:
                c = d['tracer_slices'][:]
                if c.shape != (3, meta['ny'], meta['nx']) or not np.isfinite(c).all():
                    raise ValueError(f'Invalid concentration: {file}')
                crange = [min(crange[0], float(c.min())), max(crange[1], float(c.max()))]

scale_reference_cases=[case[0]['name'] for case in cases]
if a.case_index is not None:
    cases=[cases[a.case_index]]
    a.labels=[a.labels[a.case_index]]
nrows = len(cases)*(2 if a.tracer else 1)
fig, axes = plt.subplots(nrows, 3, figsize=(15, 3.7*nrows*max(1,m["ly"]/m["lx"])), layout='constrained')
artists = []
for case_index, ((meta, files, _, _), label) in enumerate(zip(cases, a.labels)):
    fields = [('vorticity_slices', 'ζ/f₀ (asinh)', 'RdBu_r')]
    if a.tracer:
        fields.append(('tracer_slices', 'tracer C · contrast at 0.0793', 'BrBG'))
    case_artists = []
    for local_row, (key, quantity, cmap) in enumerate(fields):
        row = case_index*len(fields) + local_row
        row_artists = []
        with h5py.File(files[0]) as d:
            values = d[key][:] / meta['f'] if key == 'vorticity_slices' else d[key][:]
        for col, height in enumerate(meta['video_heights_actual_m']):
            ax = axes[row, col]
            if key == 'vorticity_slices':
                norm = AsinhNorm(linear_width=a.vorticity_width, vmin=-limits[col], vmax=limits[col], clip=False)
                ticks = [v for v in [-5,-3,-1,-.3,-.1,0,.1,.3,1,3,5] if -limits[col] <= v <= limits[col]]
            else:
                center=0.0793082; width=.001
                left=np.arcsinh((center-crange[0])/width); right=np.arcsinh((crange[1]-center)/width)
                def forward(c):
                    q=np.arcsinh((c-center)/width)
                    return np.where(c<center,q/left,q/right)
                def inverse(q):
                    return center+width*np.sinh(np.where(q<0,q*left,q*right))
                norm = FuncNorm((forward,inverse),vmin=crange[0],vmax=crange[1],clip=False)
                ticks = [0,.05,center-.002,center,center+.002,.1,1]
            im = ax.imshow(values[col], origin='lower', interpolation='nearest', cmap=cmap, norm=norm,
                           extent=(0, meta['lx']/1000, -meta['ly']/2000, meta['ly']/2000))
            ax.set(xlabel='x [km]', ylabel='y [km]',
                   title=f'{label} · {quantity}\nz = {height/1000:.2f} km above bottom')
            ax.title.set_fontsize(10)
            cb = fig.colorbar(im, ax=ax, shrink=.86)
            cb.set_ticks(ticks)
            cb.set_ticklabels([f'{v:g}' for v in ticks])
            row_artists.append(im)
        case_artists.append((key, row_artists))
    artists.append(case_artists)
title = fig.suptitle(a.title + ' | Earth day 0.000', fontsize=13)
fig.canvas.draw()
fig.set_layout_engine(None)
a.output.parent.mkdir(parents=True, exist_ok=True)
writer = FFMpegWriter(fps=a.fps, codec='libx264', extra_args=[
    '-crf', str(a.crf), '-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
with writer.saving(fig, str(a.output), dpi=160):
    for i, t in enumerate(times):
        for case_index, (meta, files, _, _) in enumerate(cases):
            with h5py.File(files[i]) as d:
                for key, ims in artists[case_index]:
                    fields = d[key][:] / meta['f'] if key == 'vorticity_slices' else d[key][:]
                    for im, field in zip(ims, fields):
                        im.set_data(field)
        title.set_text(f'{a.title} | Earth day {t/86400:.3f}')
        writer.grab_frame()
fig.savefig(a.output.with_suffix('.png'), dpi=160)
plt.close(fig)
a.output.with_suffix('.json').write_text(json.dumps({
    'scale_reference_cases':scale_reference_cases, 'labels': a.labels, 'case_names': [case[0]['name'] for case in cases],
    'title': a.title, 'layout': f'{nrows} rows, three depth columns',
    'row_order': [f'{label}: {field}' for label in a.labels for field in (['vorticity','tracer'] if a.tracer else ['vorticity'])],
    'times_seconds': times.tolist(), 'heights_m': m['video_heights_actual_m'],
    'fps': a.fps, 'frames': len(times), 'crf': a.crf,
    'vorticity_normalization': {'type':'asinh', 'linear_width':a.vorticity_width,
                               'limits_by_depth':[[-float(v),float(v)] for v in limits],
                               'shared_between_cases':True, 'fixed_in_time':True, 'units':'zeta divided by each case central f0', 'f0_by_case_s_inv':[case[0]['f'] for case in cases]},
    'tracer_normalization': {'type':'two_sided_centered_asinh','center':.0793082,'linear_width':.001,'colormap':'BrBG','limits':crange,'units':'physical concentration C','center_basis':'Rounded common background concentration from Unity tracer audit; fixed, never recalculated per frame.',
                            'shared_between_cases_and_depths':True,'fixed_in_time':True} if a.tracer else None,
    'renderer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'note':'Tracer values are unchanged absolute concentrations; only their color mapping differs from the original movies. Every saved time included. No field clipping or temporal interpolation. Values on colorbars are physical C or ζ/f₀. Rendering uses nearest-neighbor display sampling.'
}, indent=2)+'\n')
print('Rendered matched comparison:', a.output, flush=True)

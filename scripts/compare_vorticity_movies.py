"""Compare two completed cases: one case per row, three depth columns.

Each depth shares one fixed vorticity scale across both cases and every frame.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from vorticity_depth_videos import collect

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('cases', nargs=2, type=Path)
p.add_argument('--labels', nargs=2, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--fps', type=int, default=8)
a = p.parse_args()
if a.fps < 1:
    p.error('fps must be positive')
cases = [collect(path, False) for path in a.cases]
m, _, times, _ = cases[0]
for other, _, other_times, _ in cases[1:]:
    for key in ('nx', 'ny', 'nz', 'lx', 'ly', 'height', 'f', 'heat_flux',
                'cooling_model', 'cooling_time', 'video_heights_actual_m'):
        if m[key] != other[key]:
            raise ValueError(f'Incompatible comparison: {key}')
    if not np.array_equal(times, other_times):
        raise ValueError('Cases must have identical saved times')
limits = np.maximum(np.max([case[3] for case in cases], axis=0), 1e-12)
fig, axes = plt.subplots(2, 3, figsize=(15, 8), layout='constrained')
artists = []
for row, ((meta, files, _, _), label) in enumerate(zip(cases, a.labels)):
    images = []
    with h5py.File(files[0]) as data:
        fields = data['vorticity_slices'][:] / meta['f']
    for level, height in enumerate(meta['video_heights_actual_m']):
        ax = axes[row, level]
        im = ax.imshow(fields[level], origin='lower', interpolation='nearest',
                       cmap='RdBu_r', vmin=-limits[level], vmax=limits[level],
                       extent=(0, meta['lx']/1000, -meta['ly']/2000, meta['ly']/2000))
        ax.set(xlabel='x [km]', ylabel='y [km]',
               title=f'{label} | ζ/f₀\nz = {height/1000:.2f} km above bottom')
        fig.colorbar(im, ax=ax, shrink=.85)
        images.append(im)
    artists.append(images)
title = fig.suptitle('6,000-km comparison | first 10-day segments | Earth day 00.000')
fig.canvas.draw()
fig.set_layout_engine(None)
a.output.parent.mkdir(parents=True, exist_ok=True)
writer = FFMpegWriter(fps=a.fps, codec='libx264', extra_args=[
    '-crf', '26', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
with writer.saving(fig, str(a.output), dpi=100):
    for i, time in enumerate(times):
        for row, (meta, files, _, _) in enumerate(cases):
            with h5py.File(files[i]) as data:
                fields = data['vorticity_slices'][:] / meta['f']
            for im, field in zip(artists[row], fields):
                im.set_data(field)
        title.set_text(f'6,000-km comparison | first 10-day segments | Earth day {time/86400:.3f}')
        writer.grab_frame()
fig.savefig(a.output.with_suffix('.png'), dpi=100)
plt.close(fig)
a.output.with_suffix('.json').write_text(json.dumps({
    'labels': a.labels, 'case_names': [case[0]['name'] for case in cases],
    'layout': 'two case rows, three depth columns', 'times_seconds': times.tolist(),
    'heights_m': m['video_heights_actual_m'], 'fps': a.fps,
    'fixed_vorticity_limits': [[-float(v), float(v)] for v in limits],
    'note': 'Shared limits between cases at each depth, fixed over time. No spatial or temporal resampling before rendering. These older cases have no passive tracer.'
}, indent=2) + '\n')
print('Rendered synchronized comparison:', a.output, flush=True)

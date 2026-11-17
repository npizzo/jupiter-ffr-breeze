"""Unity-only media preparation for the two long 512-square runs. No 256 comparison."""
from pathlib import Path
import csv, hashlib, json, subprocess, sys, tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.animation import FFMpegWriter
import imageio_ffmpeg

src, out = map(Path, sys.argv[1:])
out.mkdir(parents=True, exist_ok=False)
cases = ['selective', 'unforced']
labels = {'selective': 'Decay + anticyclone-selective viscosity', 'unforced': 'Ordinary viscous decay · existing control'}
data, metadata = {}, {}
control_old=Path('data/viscous2d-long512-v1-1848831/unforced')
control_new=Path('data/viscous2d-fast512-v1-1872263/unforced')
for case in cases:
    frames, segments = {}, []
    paths=([(n, src/f'to_{n:03d}') for n in [80,160,240,320,400]] if case=='selective' else
           [(n,control_old/f'to_{n:03d}') for n in [40,80,120,160,200]]+[(400,control_new)])
    for end,p in paths:
        c=tomllib.loads((p/'config.toml').read_text());s=tomllib.loads((p/'status.toml').read_text())
        assert s['status']=='COMPLETED' and abs(s['t']-end)<1e-7
        assert c['grid']==[512,512] and c['viscosity']==.0001 and c['drag']==0
        assert c['f']==1 and c['beta']==0 and c['instability_coefficient']==0
        assert c['case']==('anticyclone_diffusion' if case=='selective' else 'unforced')
        assert c.get('selective_viscosity',0)==(.0004 if case=='selective' else 0)
        for row in csv.DictReader((p/'frames.csv').open()):
            it,t,path=int(row['iteration']),float(row['t']),p/row['file']
            assert abs(t-it*.0025)<1e-7
            if it in frames: assert path.read_bytes()==frames[it][1].read_bytes()
            frames[it]=(t,path)
        segments.append({'directory':str(p),'config':c,'status':s})
    assert sorted(frames)==list(range(0,160001,200))
    data[case]=[frames[it] for it in sorted(frames)];metadata[case]=segments
assert len({s['config']['initial_sha256'] for cs in metadata.values() for s in cs}) == 1
def field(case, i):
    q = np.fromfile(data[case][i][1], dtype=np.float32)
    assert q.size == 512 * 512 and np.isfinite(q).all()
    return q.reshape((512, 512), order='F')
assert np.array_equal(field('selective', 0), field('unforced', 0))
limit = max(float(np.abs(field(c, i)).max()) for c in cases for i in range(801))
assert limit > 0
norm = AsinhNorm(linear_width=.1, vmin=-limit, vmax=limit)
ff = imageio_ffmpeg.get_ffmpeg_exe()
plt.rcParams['animation.ffmpeg_path'] = ff
movies = {}
for selected, stem in [(cases, 'selective_control_512_0_400'), (['selective'], 'selective_512_0_400')]:
    fig, axs = plt.subplots(1, len(selected), figsize=(13 if len(selected) == 2 else 7.6, 6.8), squeeze=False, layout='constrained')
    images = []
    for ax, case in zip(axs.flat, selected):
        images.append(ax.imshow(field(case, 0).T, origin='lower', extent=(0, 2*np.pi, 0, 2*np.pi), cmap='RdBu_r', norm=norm, interpolation='nearest'))
        ax.set(title=labels[case], xlabel='x', ylabel='y')
    fig.colorbar(images[0], ax=list(axs.flat), shrink=.78, label='Vertical vorticity ζ · fixed shared asinh scale')
    title = fig.suptitle('')
    fig.supxlabel('512 × 512 · no forcing or drag · base viscosity = 0.0001 · one horizontal 2D plane\nPositive ζ is counterclockwise · dimensionless model time', fontsize=10)
    movie = out / f'{stem}.mp4'
    writer = FFMpegWriter(fps=12, codec='libx264', extra_args=['-b:v', '2200k', '-maxrate', '2400k', '-bufsize', '4800k', '-threads', '1', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
    with writer.saving(fig, str(movie), dpi=100):
        for i in range(801):
            for im, case in zip(images, selected):
                im.set_data(field(case, i).T)
            title.set_text(f'512² anticyclone-selective dissipation / control · t = {i * .5:.1f} / 400' if len(selected) == 2 else f'512² {labels[selected[0]].lower()} · t = {i * .5:.1f} / 400')
            writer.grab_frame()
            if i == 800:
                fig.savefig(out / f'{stem}.png', dpi=100)
    plt.close(fig)
    count, duration = imageio_ffmpeg.count_frames_and_secs(str(movie))
    assert count == 801 and abs(duration - 801 / 12) < .1
    subprocess.run([ff, '-v', 'error', '-xerror', '-i', str(movie), '-f', 'null', '-'], check=True)
    assert movie.stat().st_size < 25 * 1024**2
    movies[stem] = {'frames': count, 'duration_seconds': duration, 'full_decode_passed': True}
    print('VERIFIED', movie.name, flush=True)
report = {
    'status': 'MEDIA_VERIFIED', 'compute_location': 'Unity', 'time_range': [0, 400],
    'time_units': 'dimensionless model time, not Earth or Jupiter days',
    'frame_interval': .5, 'unique_frames_per_case': 801, 'duplicate_restart_frames_equal': True,
    'identical_initial_vorticity': True, 'segments': metadata, 'movies': movies,
    'color_scale': {'type': 'symmetric asinh', 'linear_width': .1, 'limits': [-limit, limit], 'fixed_across_cases_and_time': True, 'clipping': False},
    'notes': ['Only 512-square output; no comparison to 256-square runs.', 'No tracer or vertical layers exist in this horizontal 2D model.', 'The selective operator explicitly imposes sign asymmetry. It is not a model of convection or a centrifugal-instability threshold.', 'This is a media publication, not an inverse-cascade or convergence diagnosis.', 'Extra viscosity 0.0004 acts through Laplacian(min(zeta,0)), preserving mean circulation. Interface fluxes also affect adjacent positive vorticity. Existing unforced output is reused.'],
    'assets': {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
}
(out / 'report-private.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
public={k:v for k,v in report.items() if k not in ('segments','assets')}
(out/'report.json').write_text(json.dumps(public,indent=2,allow_nan=False)+'\n')
print('MEDIA_COMPLETE', flush=True)

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

src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
cases=['halfscale'];labels={'halfscale':'Half-wavelength instability forcing · k = 24–32'}
frames,segments={},[]
for end in [80,160,240,320,400]:
    p=src/f'to_{end:03d}'
    c=tomllib.loads((p/'config.toml').read_text());s=tomllib.loads((p/'status.toml').read_text())
    assert s['status']=='COMPLETED' and abs(s['t']-end)<1e-7
    assert c['case']=='halfscale' and c['forcing_band']==[24.,32.] and c['instability_coefficient']==.000175
    assert c['grid']==[512,512] and c['viscosity']==.0001 and c['dt']==.00125 and c['drag']==0 and c['beta']==0 and c['f']==1
    for row in csv.DictReader((p/'frames.csv').open()):
        it,t,path=int(row['iteration']),float(row['t']),p/row['file']
        assert abs(t-it*.00125)<1e-7
        if it in frames: assert path.read_bytes()==frames[it][1].read_bytes()
        frames[it]=(t,path)
    segments.append({'directory':str(p),'config':c,'status':s})
assert sorted(frames)==list(range(0,320001,400))
data={'halfscale':[frames[it] for it in sorted(frames)]};metadata={'halfscale':segments}
def field(case, i):
    q = np.fromfile(data[case][i][1], dtype=np.float32)
    assert q.size == 512 * 512 and np.isfinite(q).all()
    return q.reshape((512, 512), order='F')

limit = max(float(np.abs(field(c, i)).max()) for c in cases for i in range(801))
assert limit > 0
norm = AsinhNorm(linear_width=.1, vmin=-limit, vmax=limit)
ff = imageio_ffmpeg.get_ffmpeg_exe()
plt.rcParams['animation.ffmpeg_path'] = ff
movies = {}
for selected, stem in [(cases, 'halfscale_512_0_400')]:
    fig, axs = plt.subplots(1, len(selected), figsize=(18 if len(selected) == 3 else 7.6, 6.8), squeeze=False, layout='constrained')
    images = []
    for ax, case in zip(axs.flat, selected):
        images.append(ax.imshow(field(case, 0).T, origin='lower', extent=(0, 2*np.pi, 0, 2*np.pi), cmap='RdBu_r', norm=norm, interpolation='nearest'))
        ax.set(title=labels[case], xlabel='x', ylabel='y')
    fig.colorbar(images[0], ax=list(axs.flat), shrink=.78, label='Vertical vorticity ζ · fixed shared asinh scale')
    title = fig.suptitle('')
    fig.supxlabel('512 × 512 · ordinary viscosity only · no drag · one horizontal 2D plane\nPositive ζ is counterclockwise · dimensionless model time', fontsize=10)
    movie = out / f'{stem}.mp4'
    writer = FFMpegWriter(fps=12, codec='libx264', extra_args=['-b:v', '2200k', '-maxrate', '2400k', '-bufsize', '4800k', '-threads', '1', '-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
    with writer.saving(fig, str(movie), dpi=100):
        for i in range(801):
            for im, case in zip(images, selected):
                im.set_data(field(case, i).T)
            title.set_text(f'512² forcing strength comparison · t = {i * .5:.1f} / 400' if len(selected) == 3 else f'512² {labels[selected[0]].lower()} · t = {i * .5:.1f} / 400')
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
    'same_original_initial_velocity_file': True, 'segments': metadata, 'movies': movies,
    'color_scale': {'type': 'symmetric asinh', 'linear_width': .1, 'limits': [-limit, limit], 'fixed_across_cases_and_time': True, 'clipping': False},
    'notes': ['Only 512-square output; no comparison to 256-square runs.', 'No tracer or vertical layers exist in this horizontal 2D model.', 'Constant f cannot itself select vorticity sign in horizontal incompressible 2D dynamics.', 'This is a media publication, not an inverse-cascade or convergence diagnosis.', 'Forcing band24–32 and coefficient0.000175 match continuum net growth to band12–16 with coefficient0.0004 at viscosity0.0001. Native discrete growth differences are documented in validation. Initial turbulence remains in band12–16; it is not rescaled. This movie has no control panels.'],
    'assets': {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
}
(out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
print('MEDIA_COMPLETE', flush=True)

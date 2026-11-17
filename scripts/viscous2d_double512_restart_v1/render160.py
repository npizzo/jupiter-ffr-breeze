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

src, continuation, doubled, out = map(Path, sys.argv[1:])
out.mkdir(parents=True, exist_ok=False)
cases = ['unforced', 'forced', 'double_forced']
labels = {'forced': 'Original forcing: a = 0.0002', 'unforced': 'Unforced decay', 'double_forced': 'Double forcing: a = 0.0004'}
data, metadata = {}, {}
for case in cases:
    frames, segments = {}, []
    paths = [(160, doubled / 'to_200')] if case == 'double_forced' else [(n, src / case / f'to_{n:03d}') for n in [40, 80, 120, 160]]
    for end, p in paths:
        c = tomllib.loads((p / 'config.toml').read_text())
        s = {'status': 'PARTIAL_CFL_STOP', 'published_end_time': 160, 'last_saved_time': 161, 'stop_cfl': 0.5000724924099372} if case == 'double_forced' else tomllib.loads((p / 'status.toml').read_text())
        if case != 'double_forced':
            assert s['status'] == 'COMPLETED' and abs(s['t'] - end) < 1e-7
        assert c['grid'] == [512, 512] and c['viscosity'] == .0001 and c['drag'] == 0
        assert c['case'] == case and c['f'] == 1 and c['beta'] == 0
        assert c['instability_coefficient'] == {'forced': .0002, 'unforced': 0, 'double_forced': .0004}[case]
        for row in csv.DictReader((p / 'frames.csv').open()):
            it, t, path = int(row['iteration']), float(row['t']), p / row['file']
            assert abs(t - it * .0025) < 1e-7
            if it > 64000:
                continue
            if it in frames:
                assert path.read_bytes() == frames[it][1].read_bytes(), 'Restart frame mismatch'
            frames[it] = (t, path)
        segments.append({'directory': str(p), 'config': c, 'status': s})
    assert sorted(frames) == list(range(0, 64001, 200))
    data[case] = [frames[it] for it in sorted(frames)]
    metadata[case] = segments
assert len({s['config']['initial_sha256'] for cs in metadata.values() for s in cs}) == 1
def field(case, i):
    q = np.fromfile(data[case][i][1], dtype=np.float32)
    assert q.size == 512 * 512 and np.isfinite(q).all()
    return q.reshape((512, 512), order='F')
assert all(np.array_equal(field(c, 0), field('unforced', 0)) for c in cases)
limit = max(float(np.abs(field(c, i)).max()) for c in cases for i in range(321))
assert limit > 0
norm = AsinhNorm(linear_width=.1, vmin=-limit, vmax=limit)
ff = imageio_ffmpeg.get_ffmpeg_exe()
plt.rcParams['animation.ffmpeg_path'] = ff
movies = {}
for selected, stem in [(cases, 'forcing_strength_512_0_160'), (['double_forced'], 'double_forced_512_0_160')]:
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
        for i in range(321):
            for im, case in zip(images, selected):
                im.set_data(field(case, i).T)
            title.set_text(f'512² forcing strength comparison · t = {i * .5:.1f} / 160' if len(selected) == 3 else f'512² {labels[selected[0]].lower()} · t = {i * .5:.1f} / 160')
            writer.grab_frame()
            if i == 320:
                fig.savefig(out / f'{stem}.png', dpi=100)
    plt.close(fig)
    count, duration = imageio_ffmpeg.count_frames_and_secs(str(movie))
    assert count == 321 and abs(duration - 321 / 12) < .1
    subprocess.run([ff, '-v', 'error', '-xerror', '-i', str(movie), '-f', 'null', '-'], check=True)
    assert movie.stat().st_size < 25 * 1024**2
    movies[stem] = {'frames': count, 'duration_seconds': duration, 'full_decode_passed': True}
    print('VERIFIED', movie.name, flush=True)
report = {
    'status': 'MEDIA_VERIFIED', 'campaign_status': 'PARTIAL_THROUGH_T160', 'compute_location': 'Unity', 'time_range': [0, 160],
    'time_units': 'dimensionless model time, not Earth or Jupiter days',
    'frame_interval': .5, 'unique_frames_per_case': 321, 'duplicate_restart_frames_equal': True,
    'identical_initial_vorticity': True, 'segments': metadata, 'movies': movies,
    'color_scale': {'type': 'symmetric asinh', 'linear_width': .1, 'limits': [-limit, limit], 'fixed_across_cases_and_time': True, 'clipping': False},
    'notes': ['Only 512-square output; no comparison to 256-square runs.', 'No tracer or vertical layers exist in this horizontal 2D model.', 'Constant f cannot itself select vorticity sign in horizontal incompressible 2D dynamics.', 'This is a media publication, not an inverse-cascade or convergence diagnosis.', 'The doubled-forcing case stopped near t=161 at the CFL guard. This publication shows only t=0–160. A separate restart halves dt from t=160; no restart frames are included here. Existing controls are reused; all cases share their initial field and physical parameters except forcing strength.'],
    'assets': {p.name: {'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
}
(out / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
print('MEDIA_COMPLETE', flush=True)

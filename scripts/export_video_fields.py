"""Copy saved video planes losslessly into smaller HDF5 files for rendering.

Run against a completed segment. No spatial or temporal resampling is performed.
The original simulation output is never modified.
"""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('source', type=Path)
p.add_argument('destination', type=Path)
a = p.parse_args()
if a.destination.exists() and any(a.destination.iterdir()):
    raise ValueError('Destination must be empty')
a.destination.mkdir(parents=True, exist_ok=True)
for name in ('metadata.toml', 'diagnostics.csv'):
    shutil.copy2(a.source / name, a.destination / name)
fields = ('t', 'x', 'y', 'z', 'video_heights', 'vorticity_slices', 'tracer_slices',
          'condensate_slices', 'theta_profile', 'mean_N2_profile', 'w_rms',
          'kinetic_energy', 'cap_min', 'cap_mean')
frames = sorted(a.source.glob('frame_*.jld2'))
for file in frames:
    with h5py.File(file) as src, h5py.File(a.destination / file.name, 'w') as dst:
        for key in fields:
            if key in src:
                data = src[key][()]
                kwargs = dict(compression='gzip', compression_opts=1, shuffle=True) if data.shape else {}
                dst.create_dataset(key, data=data, **kwargs)
(a.destination / 'export.json').write_text(json.dumps({
    'source': str(a.source), 'frames': len(frames), 'resampling': False,
    'fields': fields, 'format': 'HDF5 rendering subset; not a restart file',
    'exporter_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}, indent=2) + '\n')
print(f'Exported {len(frames)} frames to {a.destination}', flush=True)

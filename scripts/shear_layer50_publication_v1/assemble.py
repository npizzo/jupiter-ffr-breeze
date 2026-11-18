"""Join saved 0–10 and 10–50 output via Unity-local symlinks; never copy raw fields locally."""
from pathlib import Path
import sys,tomllib,json,shutil
import numpy as np,h5py
old,new,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
record=[]
for cap in ('stable','neutral'):
 for shear in ('cyclonic','none','anticyclonic'):
  case=f'{cap}_{shear}';a,b=old/case,new/case
  ma,mb=[tomllib.loads((d/'metadata.toml').read_text()) for d in (a,b)]
  assert ma['status']==mb['status']=='COMPLETE'
  assert np.isclose(ma['simulated_seconds'],10*86400) and np.isclose(mb['simulated_seconds'],50*86400)
  for k in ('nx','ny','nz','lx','ly','height','f','cap_n2','shear_rossby','heat_flux','turbulence_seed','script_sha256','video_heights_actual_m'):
   assert ma[k]==mb[k],(case,k)
  dest=out/case;dest.mkdir();shutil.copyfile(b/'metadata.toml',dest/'metadata.toml')
  by_time={};duplicates=[]
  for source in (a,b):
   for p in sorted(source.glob('frame_*.jld2')):
    with h5py.File(p) as h:t=float(h['t'][()])
    if t in by_time:
     assert np.isclose(t,10*86400), 'Unexpected duplicate frame'
     with h5py.File(p) as newer,h5py.File(by_time[t]) as older:
      for key in ('vorticity_slices','tracer_slices'):
       np.testing.assert_allclose(newer[key][:],older[key][:],rtol=1e-6,atol=1e-10)
     duplicates.append(t);continue
    by_time[t]=p
  expected=np.arange(0,50*86400+1,ma['output_seconds'])
  np.testing.assert_allclose(sorted(by_time),expected,rtol=0,atol=1e-6)
  for i,t in enumerate(sorted(by_time)):(dest/f'frame_{i:05d}.jld2').symlink_to(by_time[t].resolve())
  record.append(dict(case=case,frames=len(expected),duplicate_boundary_times=duplicates,source_segments=[str(a),str(b)]))
(out/'assembly.json').write_text(json.dumps(record,indent=2)+'\n')
print('FIFTY_DAY_TIMELINES_VERIFIED',flush=True)

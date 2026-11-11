"""Render missing earlier-campaign movies on Unity; never change simulation output."""
from pathlib import Path
import json, tomllib, subprocess, sys, shutil, hashlib, os
import imageio_ffmpeg
import h5py
root=Path.cwd(); out=root/'rendered/earlier-campaigns-20261111'; out.mkdir(exist_ok=True)
items=[
('polar_unheated_2days','data/convection2026-v1/convection2026_f0/segment_1','85 N unheated pilot | completed days 0-2',False),
('polar_pilot_f0','data/convection2026-preflight-v1/pilot512_f0','85 N validation | 0 W/m² | six hours',False),
('polar_pilot_f54','data/convection2026-preflight-v1/pilot512_f54','85 N validation | 54 W/m² | six hours',False),
('polar_pilot_f162','data/convection2026-preflight-v1/pilot512_f162','85 N validation | 162 W/m² | six hours',False),
('shear_fine_partial','data/shear20-v2/shear20_fine/segment_1','Fine cyclonic-shear run | saved partial output',True),
('moist_rest','data/moist-validation-v3/moist_rest','Moist validation | resting atmosphere | one hour',False),
('shear_grid1024','data/shear20-preflight-v4/grid_probe_1024','Shear validation | 1024² grid | two hours',False),
('shear_grid1536','data/shear20-preflight-v4/grid_probe_1536','Shear validation | 1536² grid | two hours',False),
('shear_balanced_rest','data/shear20-preflight-v4/balanced_rest','Shear validation | balanced rest | one day',False),
('shear_uniform_tracer','data/shear20-preflight-v4/passivity_uniform','Shear validation | uniform passive tracer | one hour',False),
('large_grid_probe','data/large40-preflight-v2/large_grid_probe','Large-domain validation | two hours',False),
]
report={'compute_location':'Unity','raw_fields_downloaded':False,'status':'RENDERING','movies':[]}
for name,source,label,preview in items:
 src=root/source; meta=tomllib.loads((src/'metadata.toml').read_text())
 view=out/'views'/name; view.mkdir(parents=True,exist_ok=True)
 shutil.copyfile(src/'metadata.toml',view/'metadata.toml')
 frames=sorted(src.glob('frame_*.jld2'))
 # Omit an actively written last frame; a later completion can supersede this preview.
 if preview and meta['status']!='COMPLETE': frames=frames[:-1]
 for frame in frames:
  link=view/frame.name
  if not link.exists():link.symlink_to(frame)
 target=out/(name+'.mp4')
 args=[str(root/'scripts/flow_tracer_video.py'),str(view),'--output',str(target),'--label',label,'--vorticity-linear-width','0.1','--crf','30']
 if preview:args+=['--preview']
 with h5py.File(frames[0]) as d:
  if 'tracer_slices' not in d:args+=['--vorticity-only']
  elif name=='moist_rest':args+=['--cloud']
 command=[sys.executable,'-c','import matplotlib, imageio_ffmpeg, runpy, sys; matplotlib.rcParams["animation.ffmpeg_path"]=imageio_ffmpeg.get_ffmpeg_exe(); sys.argv=sys.argv[1:]; runpy.run_path(sys.argv[0],run_name="__main__")']+args
 subprocess.run(command,check=True,env={**os.environ,'PYTHONPATH':str(root/'scripts')})
 m=json.loads(target.with_suffix('.json').read_text()); count,duration=imageio_ffmpeg.count_frames_and_secs(str(target))
 assert count==len(m['times_seconds']) and target.stat().st_size<24000000
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(target),'-f','null','-'],check=True)
 shutil.copyfile(src/'metadata.toml',out/(name+'_metadata.toml'))
 report['movies'].append({'name':name,'source':source,'preview':preview,'start_day':m['times_seconds'][0]/86400,'end_day':m['times_seconds'][-1]/86400,'frames':count,'duration_seconds':duration,'full_decode_passed':True,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()})
 (out/'publication.json').write_text(json.dumps(report,indent=2)+'\n'); print('VERIFIED',report['movies'][-1],flush=True)
report['status']='COMPLETE'; (out/'publication.json').write_text(json.dumps(report,indent=2)+'\n')

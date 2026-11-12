"""Completed heated 35N beta above heated 55N beta, four rows by three depths; Unity only."""
from pathlib import Path
import sys,json,tomllib,hashlib,subprocess,shutil,os
import imageio_ffmpeg,h5py,numpy as np
root=Path.cwd()
out=root/'rendered/latitude100-contrast';out.mkdir(exist_ok=False)
cases=[root/'data/cascade35-v1/cascade35_heated',root/'data/cascadebeta100-v1/cascadebeta100_heated']
report={'status':'PREPARING','compute_location':'Unity','raw_fields_downloaded':False,'cases':{},'media':{}}
reference_config=tomllib.loads((root/'configs/cascadebeta100_heated.toml').read_text())
for label,case in zip(['latitude35','latitude55'],cases):
 records=[]
 for seg in range(1,11):
  p=case/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text())
  assert m['status']=='COMPLETE' and m['stop_seconds']==seg*10*86400
  for key,value in reference_config.items():
   if key not in ['name','stop_seconds','latitude','beta','turbulence_rossby']:assert m[key]==value,(label,seg,key)
  assert m['rotation_model']=='beta-plane'
  assert m['latitude']==(35 if label=='latitude35' else 55)
  assert np.isclose(m['turbulence_rossby']*np.sin(np.deg2rad(m['latitude'])),reference_config['turbulence_rossby']*np.sin(np.deg2rad(reference_config['latitude'])),rtol=1e-12)
  assert m['beta']==(4.02976240665303e-12 if label=='latitude35' else 2.8216700143325537e-12)
  assert json.loads((p/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
  shutil.copyfile(p/'metadata.toml',out/f'{label}_segment{seg}_metadata.toml')
  records.append({'segment':seg,'source':str(p.relative_to(root)),'frames':len(list(p.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
 assert sum(v['frames'] for v in records)==801
 report['cases'][label]=records

# Confirm the actual initial saved vorticity and tracer match, in addition to parameters.
initial=[]
for case in cases:
 with h5py.File(sorted((case/'segment_1').glob('frame_*.jld2'))[0]) as q:
  initial.append({k:q[k][:] for k in ['u','v','w','theta_profile','vorticity_slices','tracer_slices']})
for k in initial[0]:assert np.allclose(initial[0][k],initial[1][k],rtol=1e-12,atol=1e-13),'Initial fields differ: '+k
target=out/'latitude35_over55_0_100days.mp4'
command=[sys.executable,str(Path(__file__).with_name('render_movies.py')),*map(str,cases),'--labels','35°N heated beta-plane','55°N heated beta-plane','--title','Jupiter | 35°N above 55°N | matched initial eddies and bottom heating | days 0–100','--tracer','--fps','12','--crf','29','--output',str(target)]
subprocess.run(command,check=True)
m=json.loads(target.with_suffix('.json').read_text())
assert m['times_seconds']==[i*10800 for i in range(801)] and m['layout']=='4 rows, three depth columns'
assert m['row_order']==['35°N heated beta-plane: vorticity','35°N heated beta-plane: tracer','55°N heated beta-plane: vorticity','55°N heated beta-plane: tracer']
original=target.with_name(target.stem+'_original.mp4');shutil.copyfile(target,original)
for crf in [31,33,35,37,39]:
 if target.stat().st_size<23_000_000:break
 temp=target.with_name(target.stem+'_web.mp4')
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 temp.replace(target);m['web_reencode_crf']=crf
assert target.stat().st_size<24_000_000
count,duration=imageio_ffmpeg.count_frames_and_secs(str(target));assert count==801 and abs(duration-801/12)<.03
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(target),'-f','null','-'],check=True)
qa=[]
for day in [0,50,100]:
 q=out/f'qa_day{day}.jpg'
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-ss',str(day*8/12),'-i',str(target),'-frames:v','1','-q:v','2',str(q)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);qa.append(q.name)
m.update(full_decode_passed=True,compute_location='Unity',initial_saved_fields_match=True,initial_match_tolerance="rtol=1e-12, atol=1e-13",slurm_job_id=os.environ.get('SLURM_JOB_ID'),sha256=hashlib.sha256(target.read_bytes()).hexdigest())
target.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
report.update(status='COMPLETE',initial_saved_fields_match=True,initial_match_tolerance="rtol=1e-12, atol=1e-13",slurm_job_id=os.environ.get('SLURM_JOB_ID'))
report['media'][target.name]={'frames':count,'duration_seconds':duration,'bytes':target.stat().st_size,'sha256':m['sha256'],'full_decode_passed':True}
(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
names=[target.name,target.with_suffix('.png').name,target.with_suffix('.json').name,'provenance.json']+qa+[p.name for p in out.glob('*_metadata.toml')]
assets={}
for n in names:
 raw=(out/n).read_bytes();assets[n]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'git_blob_sha1':hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','full_decode_passed':True,'end_day':100,'frames':count,'duration_seconds':duration,'compute_location':'Unity','assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

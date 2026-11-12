"""Completed 55N beta over heated f-plane, four rows by three depths; Unity only."""
from pathlib import Path
import sys,json,tomllib,hashlib,subprocess,shutil,os
import imageio_ffmpeg,h5py,numpy as np
root=Path.cwd()
out=root/'rendered/beta55-f100-contrast';out.mkdir(exist_ok=False)
cases=[root/'data/cascadebeta100-v1/cascadebeta100_heated',root/'data/cascade100-v1/cascade100_f_heated']
report={'status':'PREPARING','compute_location':'Unity','raw_fields_downloaded':False,'cases':{},'media':{}}
reference_config=tomllib.loads((root/'configs/cascade100_f_heated.toml').read_text())
for label,case in zip(['beta_plane','f_plane'],cases):
 records=[]
 for seg in range(1,11):
  p=case/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text())
  assert m['status']=='COMPLETE' and m['stop_seconds']==seg*10*86400
  for key,value in reference_config.items():
   if key not in ['name','stop_seconds','rotation_model','beta']:assert m[key]==value,(label,seg,key)
  assert m['rotation_model']==('f-plane' if label=='f_plane' else 'beta-plane')
  assert m['beta']==(0 if label=='f_plane' else 2.8216700143325537e-12)
  assert json.loads((p/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
  shutil.copyfile(p/'metadata.toml',out/f'{label}_segment{seg}_metadata.toml')
  records.append({'segment':seg,'source':str(p.relative_to(root)),'frames':len(list(p.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
 assert sum(v['frames'] for v in records)==801
 report['cases'][label]=records

# Confirm the actual initial saved vorticity and tracer match, in addition to parameters.
initial=[]
for case in cases:
 with h5py.File(sorted((case/'segment_1').glob('frame_*.jld2'))[0]) as q:
  initial.append({k:q[k][:] for k in ['vorticity_slices','tracer_slices']})
for k in initial[0]:assert np.array_equal(initial[0][k],initial[1][k]),'Initial fields differ: '+k
target=out/'beta_over_f_0_100days.mp4'
command=[sys.executable,str(Path(__file__).with_name('render_movies.py')),*map(str,cases),'--labels','55°N heated beta-plane','55°N heated f-plane','--title','Jupiter 55°N | beta above f-plane | identical initial eddies and heating | days 0–100','--tracer','--fps','12','--crf','29','--output',str(target)]
subprocess.run(command,check=True)
m=json.loads(target.with_suffix('.json').read_text())
assert m['times_seconds']==[i*10800 for i in range(801)] and m['layout']=='4 rows, three depth columns'
assert m['row_order']==['55°N heated beta-plane: vorticity','55°N heated beta-plane: tracer','55°N heated f-plane: vorticity','55°N heated f-plane: tracer']
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
m.update(full_decode_passed=True,compute_location='Unity',initial_saved_fields_identical=True,slurm_job_id=os.environ.get('SLURM_JOB_ID'),sha256=hashlib.sha256(target.read_bytes()).hexdigest())
target.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
report.update(status='COMPLETE',initial_saved_fields_identical=True,slurm_job_id=os.environ.get('SLURM_JOB_ID'))
report['media'][target.name]={'frames':count,'duration_seconds':duration,'bytes':target.stat().st_size,'sha256':m['sha256'],'full_decode_passed':True}
(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
names=[target.name,target.with_suffix('.png').name,target.with_suffix('.json').name,'provenance.json']+qa+[p.name for p in out.glob('*_metadata.toml')]
assets={}
for n in names:
 raw=(out/n).read_bytes();assets[n]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest(),'git_blob_sha1':hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','full_decode_passed':True,'end_day':100,'frames':count,'duration_seconds':duration,'compute_location':'Unity','assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

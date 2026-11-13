"""Publish all saved days 0–200 of the heated 35N run; Unity only."""
from pathlib import Path
import json,tomllib,hashlib,subprocess,sys,shutil,os
import numpy as np,imageio_ffmpeg
root=Path.cwd();out=root/'rendered/beta35-full200';out.mkdir(exist_ok=False)
case=out/'input_case';case.mkdir()
cfg=tomllib.loads((root/'configs/cascade35_heated.toml').read_text())
records=[]
for seg in range(1,21):
 source=root/('data/cascade35-v1' if seg<=10 else 'data/cascade35-300-v1')/'cascade35_heated'/f'segment_{seg}'
 m=tomllib.loads((source/'metadata.toml').read_text())
 assert m['status']=='COMPLETE' and m['simulated_seconds']==seg*10*86400
 for k,v in cfg.items():
  if k!='stop_seconds':assert m[k]==v,(seg,k)
 assert json.loads((source/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
 (case/f'segment_{seg}').symlink_to(source,target_is_directory=True)
 shutil.copyfile(source/'metadata.toml',out/f'segment{seg}_metadata.toml')
 records.append({'segment':seg,'source':str(source.relative_to(root)),'frames':len(list(source.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
assert sum(x['frames'] for x in records)==1601
movie=out/'latitude35_0_200days.mp4'
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),str(case),'--labels','35°N heated beta-plane','--title','Jupiter 35°N | heated beta-plane | days 0–200','--tracer','--fps','12','--crf','31','--output',str(movie)],check=True)
m=json.loads(movie.with_suffix('.json').read_text())
assert m['times_seconds']==[i*10800 for i in range(1601)]
assert m['layout']=='2 rows, three depth columns'
assert m['row_order']==['35°N heated beta-plane: vorticity','35°N heated beta-plane: tracer']
original=movie.with_name(movie.stem+'_original.mp4');shutil.copyfile(movie,original)
for crf in [33,35,37,39]:
 if movie.stat().st_size<23_000_000:break
 temp=movie.with_name(movie.stem+'_web.mp4')
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 temp.replace(movie);m['web_reencode_crf']=crf
assert movie.stat().st_size<24_000_000
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==1601 and abs(duration-1601/12)<.04
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
qa=[]
for day in [0,100,200]:
 q=out/f'qa_day{day}.jpg'
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-ss',str(day*8/12),'-i',str(movie),'-frames:v','1','-q:v','2',str(q)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);qa.append(q.name)
m.update(full_decode_passed=True,compute_location='Unity',slurm_job_id=os.environ.get('SLURM_JOB_ID'),sha256=hashlib.sha256(movie.read_bytes()).hexdigest())
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
report={'status':'COMPLETE','compute_location':'Unity','source_segments':records,'start_day':0,'end_day':200,'original_day100_campaign':'data/cascade35-v1','continuation_campaign':'data/cascade35-300-v1','frames':count,'duration_seconds':duration,'full_decode_passed':True,'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'raw_fields_downloaded':False,'renderer_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
names=[movie.name,movie.with_suffix('.png').name,movie.with_suffix('.json').name,'provenance.json']+qa+[p.name for p in out.glob('*_metadata.toml')]
assets={}
for name in names:
 b=(out/name).read_bytes();assets[name]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'git_blob_sha1':hashlib.sha1(('blob '+str(len(b))+'\0').encode()+b).hexdigest()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','full_decode_passed':True,'frames':count,'duration_seconds':duration,'end_day':200,'assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

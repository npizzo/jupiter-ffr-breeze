"""Render completed Ferrel width campaign assets on Unity, without modifying fields."""
from pathlib import Path
import json,tomllib,hashlib,subprocess,sys,shutil,os
import imageio_ffmpeg
mode=sys.argv[1];root=Path(sys.argv[2]);out=Path(sys.argv[3]);assert mode in ['diagnostic','wide'];out.mkdir(exist_ok=False,parents=True)
case=out/'input_case';case.mkdir();records=[]
start,end=(200,210) if mode=='diagnostic' else (0,200)
segments=[root/mode/'segment_1'] if mode=='diagnostic' else [root/mode/f'segment_{i}' for i in range(1,21)]
for i,source in enumerate(segments,1):
 m=tomllib.loads((source/'metadata.toml').read_text());assert m['status']=='COMPLETE'
 assert m['simulated_seconds']==(210 if mode=='diagnostic' else i*10)*86400
 records.append({'source':str(source),'metadata_sha256':hashlib.sha256((source/'metadata.toml').read_bytes()).hexdigest()})
 if mode=='diagnostic':
  target=case/'segment_1';target.mkdir();shutil.copyfile(source/'metadata.toml',target/'metadata.toml')
  for p in source.glob('frame_*.jld2'):(target/p.name).symlink_to(p.resolve())
  first=Path('data/cascade35-300-v1/cascade35_heated/segment_20/frame_01600.jld2')
  assert first.is_file()
  if not (target/first.name).exists():(target/first.name).symlink_to(first.resolve())
 else:(case/f'segment_{i}').symlink_to(source.resolve(),target_is_directory=True)
label='35°N diagnostic restart' if mode=='diagnostic' else '35°N wider channel'
movie=out/f'{mode}_{start}_{end}days.mp4'
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),str(case),'--labels',label,'--title',f'{label} | Earth days {start}–{end}','--tracer','--fps','12','--crf','33','--output',str(movie)],check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(start*8,end*8+1)]
original=out/'publication_original.mp4';shutil.copyfile(movie,original)
for crf in [35,37,39,41]:
 if movie.stat().st_size<23_000_000:break
 tmp=out/'web.mp4';subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(tmp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);tmp.replace(movie);m['web_reencode_crf']=crf
assert movie.stat().st_size<24_000_000
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==(end-start)*8+1
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
m.update(full_decode_passed=True,compute_location='Unity',slurm_job=os.environ.get('SLURM_JOB_ID'));movie.with_suffix('.json').write_text(json.dumps(m,indent=2))
assets={}
for p in [movie,movie.with_suffix('.json'),movie.with_suffix('.png')]:
 b=p.read_bytes();assets[p.name]={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','frames':frames,'duration':duration,'start_day':start,'end_day':end,'sources':records,'assets':assets,'full_decode_passed':True},indent=2)+'\n')
print('READY_FOR_GALLERY_PUBLICATION',out,flush=True)

"""Render one of three contrast movies on Unity; retain original media unchanged."""
from pathlib import Path
import sys,os,json,tomllib,hashlib,subprocess,shutil
import imageio_ffmpeg
root=Path.cwd();out=root/'rendered/tracer100-contrast';out.mkdir(exist_ok=True)
index=int(os.environ['SLURM_ARRAY_TASK_ID']);name=['comparison','heated','unheated'][index]
audit=json.loads((out/'audit.json').read_text());assert audit['status']=='COMPLETE'
cases=[root/'data/cascade100-v1/cascade100_f_heated',root/'data/cascade100-v1/cascade100_f_control']
for case in cases:
 for seg in range(1,11):
  m=tomllib.loads((case/f'segment_{seg}'/'metadata.toml').read_text());assert m['status']=='COMPLETE' and m['stop_seconds']==seg*10*86400
movie=out/f'{name}_contrast_0_100days.mp4'
cmd=[sys.executable,str(Path(__file__).with_name('render_movies.py')),*map(str,cases),'--labels','Heating ON','Heating OFF','--title','Days 0–100 | tracer contrast centered at C=0.0793 | fixed scales, unchanged fields','--tracer','--fps','12','--crf','29','--output',str(movie)]
if index:cmd+=['--case-index',str(index-1)]
subprocess.run(cmd,check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(801)]
reference=root/'rendered/cascade100-full100-gallery'/('comparison_0_100days.json')
old=json.loads(reference.read_text());assert m['vorticity_normalization']==old['vorticity_normalization']
archive=movie.with_name(movie.stem+'_original.mp4');shutil.copyfile(movie,archive)
for crf in [31,33,35,37,39]:
 if movie.stat().st_size<=10_000_000:break
 temp=movie.with_name(movie.stem+'_web.mp4')
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(archive),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 temp.replace(movie);m['web_reencode_crf']=crf
assert movie.stat().st_size<24_000_000
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==801 and abs(duration-801/12)<.03
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
m.update(compute_location='Unity',full_decode_passed=True,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),bytes=movie.stat().st_size,unchanged_vorticity_normalization_verified=True)
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
(out/f'PART{index}_PASS.json').write_text(json.dumps({'status':'COMPLETE','file':movie.name,'poster':movie.with_suffix('.png').name,'frames':count,'duration_seconds':duration,'sha256':m['sha256'],'bytes':m['bytes'],'full_decode_passed':True},indent=2)+'\n')
print('VERIFIED',movie.name,m['bytes'],flush=True)

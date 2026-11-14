"""Render the completed 300-day run with the established tracer contrast, on Unity."""
from pathlib import Path
import subprocess,sys,json,hashlib,os
import imageio_ffmpeg
out=Path('rendered/cascade35-full300'); movie=out/'heated35_0_300days_contrast.mp4'
assert not movie.exists()
subprocess.run([sys.executable,'scripts/ferrel35_width_gallery/render_movies.py','rendered/cascade35-full300-view','--labels','35°N heated beta-plane','--title','35°N heated beta-plane | Earth days 0–300 | tracer contrast','--tracer','--fps','12','--crf','33','--output',str(movie)],check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(2401)]
original=movie.with_name('heated35_0_300days_contrast_original.mp4');movie.rename(original)
for crf in [33,35,37,39]:
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(movie)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 if movie.stat().st_size<24000000:break
assert movie.stat().st_size<24000000
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==2401 and abs(duration-2401/12)<.03
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
m.update(full_decode_passed=True,compute_location='Unity',slurm_job=os.environ.get('SLURM_JOB_ID'),web_reencode_crf=crf)
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [movie,movie.with_suffix('.png'),movie.with_suffix('.json')]}
(out/'contrast_PASS.json').write_text(json.dumps({'status':'COMPLETE','start_day':0,'end_day':300,'frames':frames,'duration_seconds':duration,'compute_location':'Unity','full_decode_passed':True,'assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',movie,flush=True)

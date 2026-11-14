"""Short diagnostic movie, rendered and fully decoded on Unity."""
from pathlib import Path
import sys,subprocess,shutil,json,hashlib,os
import imageio_ffmpeg
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
case=out/'input_case';segment=case/'segment_1';segment.mkdir(parents=True,exist_ok=False)
shutil.copyfile(root/'metadata.toml',segment/'metadata.toml')
for p in root.glob('frame_*.jld2'):(segment/p.name).symlink_to(p.resolve())
p=Path('data/cascade35-300-v1/cascade35_heated/segment_30/frame_02400.jld2');(segment/p.name).symlink_to(p.resolve())
movie=out/'diagnostic_300_301days.mp4'
subprocess.run([sys.executable,'scripts/ferrel35_width_gallery/render_movies.py',str(case),'--labels','35°N native-face diagnostic','--title','35°N heated beta-plane | Earth days 300–301 | diagnostic restart','--tracer','--fps','1','--crf','28','--output',str(movie)],check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(2400,2409)]
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==9
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
m.update(full_decode_passed=True,compute_location='Unity',slurm_job=os.environ.get('SLURM_JOB_ID'));movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.glob('*') if p.is_file()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','start_day':300,'end_day':301,'frames':frames,'duration_seconds':duration,'full_decode_passed':True,'compute_location':'Unity','assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

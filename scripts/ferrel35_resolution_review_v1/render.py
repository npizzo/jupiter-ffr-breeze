"""Render matched coarse/fine movies on Unity and verify the entire decode."""
from pathlib import Path
import sys,subprocess,json,hashlib,os
import imageio_ffmpeg
root=Path(sys.argv[1]);out=Path(sys.argv[2]);assert (out/'resolution_review.json').exists()
cases=[]
for label in ['coarse','fine']:
 case=out/(label+'_view');case.mkdir(exist_ok=False)
 (case/'segment_1').symlink_to((root/label).resolve(),target_is_directory=True);cases.append(str(case))
movie=out/'coarse_over_fine_300_301days.mp4'
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),*cases,'--labels','256² × 64 | Δt ≤ 15 s','512² × 64 | Δt ≤ 15 s','--title','35°N heated beta-plane | horizontal resolution test | days 300–301','--tracer','--fps','1','--crf','28','--output',str(movie)],check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(2400,2409)]
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert frames==9 and duration==9
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
assert movie.stat().st_size<24_000_000
m.update(full_decode_passed=True,compute_location='Unity',slurm_job=os.environ.get('SLURM_JOB_ID'),spatial_comparison='Different native horizontal grids; same physical extent and levels; no field interpolation in the movie')
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','start_day':300,'end_day':301,'frames':frames,'duration_seconds':duration,'full_decode_passed':True,'assets':assets},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

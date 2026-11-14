from pathlib import Path
import sys,subprocess,json,hashlib,os
import imageio_ffmpeg
from diagnostics import analyze
root=Path(sys.argv[1]);name=sys.argv[2];out=Path(sys.argv[3])/name
analyze(root/name,out)
subprocess.run([sys.executable,str(Path(__file__).with_name('plot.py')),str(out)],check=True)
movie=out/(name+'.mp4')
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),str(root/name),'--labels',name,'--title','Jupiter f-plane sign-selection tests | '+name,'--tracer','--fps','4','--crf','28','--output',str(movie)],check=True)
frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie))
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
assert movie.stat().st_size<24_000_000
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.suffix in ['.json','.png','.mp4']}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','frames':frames,'duration_seconds':duration,'full_decode_passed':True,'slurm_job':os.environ.get('SLURM_JOB_ID'),'assets':assets},indent=2)+'\n')
print('POINT4_PUBLICATION_ASSETS_READY',name,flush=True)

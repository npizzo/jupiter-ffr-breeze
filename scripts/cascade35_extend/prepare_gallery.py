"""Join original and continued output without changing model data; render on Unity."""
from pathlib import Path
import sys,json,tomllib,hashlib,shutil,subprocess
import imageio_ffmpeg
root=Path.cwd();extension,original=[Path(p).resolve() for p in sys.argv[1:3]]
out=root/'rendered/cascade35-full300';out.mkdir(exist_ok=True)
view=root/'rendered/cascade35-full300-view';view.mkdir(exist_ok=True)
config=tomllib.loads((root/'configs/cascade35_heated.toml').read_text())
records=[]
for seg in range(1,31):
 p=(original if seg<=10 else extension)/'cascade35_heated'/f'segment_{seg}'
 m=tomllib.loads((p/'metadata.toml').read_text())
 assert m['status']=='COMPLETE' and m['simulated_seconds']==seg*10*86400
 assert m['segment_start_seconds']==(seg-1)*10*86400
 for k,v in config.items():
  if k!='stop_seconds':assert m[k]==v,(seg,k)
 assert json.loads((p/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
 if seg>10:assert tomllib.loads((p/'continuation.toml').read_text())['status']=='PASS'
 link=view/f'segment_{seg}'
 if not link.exists():link.symlink_to(p,target_is_directory=True)
 assert link.resolve()==p.resolve()
 shutil.copyfile(p/'metadata.toml',out/f'segment{seg}_metadata.toml')
 records.append({'segment':seg,'source':str(p.relative_to(root)),'frames':len(list(p.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
assert sum(r['frames'] for r in records)==2401
movie=out/'heated35_0_300days.mp4'
subprocess.run([sys.executable,str(Path(__file__).with_name('render_movie.py')),str(view),'--output',str(movie),'--label','35°N heated beta-plane | full checkpoint continuation | days 0–300','--vorticity-linear-width','0.1','--tracer-linear-width','0.01','--fps','12','--crf','29'],check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(2401)]
archive=movie.with_name(movie.stem+'_original.mp4');shutil.copyfile(movie,archive)
for crf in [31,33,35,37,39]:
 if movie.stat().st_size<=12_000_000:break
 temp=movie.with_name(movie.stem+'_web.mp4')
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(archive),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 temp.replace(movie);m['web_reencode_crf']=crf
assert movie.stat().st_size<24_000_000
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==2401 and abs(duration-2401/12)<.03
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
m.update(compute_location='Unity',frames=count,full_decode_passed=True,final_mp4_sha256=hashlib.sha256(movie.read_bytes()).hexdigest())
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
subprocess.run([sys.executable,str(Path(__file__).with_name('jet_review.py')),str(view),str(out)],check=True)
(out/'provenance.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','raw_fields_downloaded':False,'start_day':0,'end_day':300,'segments':records,'media':{movie.name:{'frames':count,'duration_seconds':duration,'sha256':m['final_mp4_sha256'],'bytes':movie.stat().st_size,'full_decode_passed':True}},'next_decision':'Review day300-review.json before any extension to day 500.'},indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

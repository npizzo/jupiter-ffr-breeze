"""Create watchable, verified archive assets from Unity renders; preserve every frame."""
from pathlib import Path
import json, shutil, subprocess, hashlib, os, sys, tomllib
import imageio_ffmpeg
root=Path.cwd(); src=root/'rendered/earlier-campaigns-20261111'; out=root/'rendered/earlier-campaigns-20261111-web';out.mkdir(exist_ok=True)
r=json.loads((src/'publication.json').read_text()); assert r['status']=='COMPLETE'
for f in src.iterdir():
 if f.is_file():shutil.copyfile(f,out/f.name)
for mode,label in [('heated','Heating ON'),('control','Heating OFF')]:
 name='cascade_day1_'+mode; case=root/'data/cascade100-preflight-v1'/('day1_'+mode);target=out/(name+'.mp4')
 command=[sys.executable,'-c','import matplotlib,imageio_ffmpeg,runpy,sys;matplotlib.rcParams["animation.ffmpeg_path"]=imageio_ffmpeg.get_ffmpeg_exe();sys.argv=sys.argv[1:];runpy.run_path(sys.argv[0],run_name="__main__")',str(root/'scripts/flow_tracer_video.py'),str(case),'--output',str(target),'--label','100-day campaign validation | '+label+' | first day','--vorticity-linear-width','0.1','--fps','1','--crf','30']
 subprocess.run(command,check=True,env={**os.environ,'PYTHONPATH':str(root/'scripts')});shutil.copyfile(case/'metadata.toml',out/(name+'_metadata.toml'))
 r['movies'].append({'name':name,'source':str(case.relative_to(root)),'preview':False,'start_day':0,'end_day':1})
for item in r['movies']:
 target=out/(item['name']+'.mp4'); meta=json.loads(target.with_suffix('.json').read_text());count=len(meta['times_seconds']); fps=max(1,min(6,count//7));oldfps=meta['fps']
 if oldfps!=fps:
  temp=target.with_name(target.stem+'_retimed.mp4')
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-v','error','-i',str(target),'-vf',f'setpts={oldfps/fps}*PTS','-r',str(fps),'-c:v','libx264','-crf','30','-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True)
  n,duration=imageio_ffmpeg.count_frames_and_secs(str(temp));assert n==count
  temp.replace(target);meta['original_render_fps']=oldfps;meta['fps']=fps;meta['playback_retiming']='Same saved frames in order; no temporal interpolation.'
 n,duration=imageio_ffmpeg.count_frames_and_secs(str(target));assert n==count
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(target),'-f','null','-'],check=True)
 item.update(frames=n,duration_seconds=duration,full_decode_passed=True,sha256=hashlib.sha256(target.read_bytes()).hexdigest())
 meta['final_mp4_sha256']=item['sha256'];target.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
r['source_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(out/'publication.json').write_text(json.dumps(r,indent=2)+'\n');print('READY',json.dumps(r),flush=True)

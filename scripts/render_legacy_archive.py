"""Render authorized single-height historical output on allocated Unity CPUs."""
from pathlib import Path
import json,os,sys,subprocess,shutil,hashlib,tomllib
import imageio_ffmpeg
root=Path.cwd();out=root/'rendered/legacy-campaigns-20261111';out.mkdir(exist_ok=True)
inv=json.loads((root/'rendered/gallery-coverage-20261111/inventory.json').read_text())['cases']
cases=[c for c in inv if c['classification']=='legacy_output_without_three_saved_vorticity_planes' and c['frame_count']>=10]
index=int(os.environ['SLURM_ARRAY_TASK_ID']);parts=4
for i,c in enumerate(cases):
 if i%parts!=index:continue
 source=root/c['source']; name=source.parent.name+'__'+source.name; target=out/(name+'.mp4');meta=tomllib.loads((source/'metadata.toml').read_text())
 label=source.name.replace('controlled48','controlled').replace('convection120','5-day').replace('convection48','convection').replace('turbulence48','turbulence').replace('neutral48','neutral').replace('strong48','strong').replace('regional20','regional').replace('_',' ')
 fps=max(1,min(12,c['frame_count']//8))
 args=[str(root/'scripts/legacy_vorticity_video.py'),str(source),'--output',str(target),'--fps',str(fps),'--labels',label]
 if meta['status']!='COMPLETE':args+=['--allow-sparse']
 command=[sys.executable,'-c','import matplotlib,imageio_ffmpeg,runpy,sys;matplotlib.rcParams["animation.ffmpeg_path"]=imageio_ffmpeg.get_ffmpeg_exe();sys.argv=sys.argv[1:];runpy.run_path(sys.argv[0],run_name="__main__")']+args
 subprocess.run(command,check=True)
 count,duration=imageio_ffmpeg.count_frames_and_secs(str(target));m=json.loads(target.with_suffix('.json').read_text());assert count==m['frames'] and target.stat().st_size<24000000
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(target),'-f','null','-'],check=True)
 m.update(source=c['source'],compute_location='Unity',full_decode_passed=True,duration_seconds=duration,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),renderer_sha256=hashlib.sha256((root/'scripts/legacy_vorticity_video.py').read_bytes()).hexdigest());target.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n');shutil.copyfile(source/'metadata.toml',out/(name+'_metadata.toml'))
 print('VERIFIED',name,count,flush=True)
(out/f'PART{index}_PASS.json').write_text(json.dumps({'status':'COMPLETE','part':index,'cases':[c['source'] for i,c in enumerate(cases) if i%parts==index]},indent=2)+'\n')

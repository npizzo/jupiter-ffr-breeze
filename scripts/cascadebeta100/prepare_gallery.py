"""Render and validate the matched heated beta/f-plane comparison on Unity."""
from pathlib import Path
import sys,json,tomllib,hashlib,subprocess,shutil
import imageio_ffmpeg
root=Path.cwd(); source=Path(sys.argv[1]).resolve()
out=root/'rendered/cascadebeta100-full100';out.mkdir(exist_ok=True)
cases=[root/'data/cascade100-v1/cascade100_f_heated',source/'cascadebeta100_heated']
report={'status':'PREPARING','compute_location':'Unity','raw_fields_downloaded':False,'cases':{},'media':{}}
reference_config=tomllib.loads((root/'configs/cascade100_f_heated.toml').read_text())
for label,case in zip(['f_plane','beta_plane'],cases):
 records=[]
 for seg in range(1,11):
  p=case/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text())
  assert m['status']=='COMPLETE' and m['stop_seconds']==seg*10*86400
  for key,value in reference_config.items():
   if key not in ['name','stop_seconds','rotation_model','beta']:assert m[key]==value,(label,seg,key)
  assert m['rotation_model']==('f-plane' if label=='f_plane' else 'beta-plane')
  assert m['beta']==(0 if label=='f_plane' else 2.8216700143325537e-12)
  assert json.loads((p/'scale_transfer/provenance.json').read_text())['status']=='COMPLETE'
  shutil.copyfile(p/'metadata.toml',out/f'{label}_segment{seg}_metadata.toml')
  records.append({'segment':seg,'source':str(p.relative_to(root)),'frames':len(list(p.glob('frame_*.jld2'))),'driver_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
 assert sum(v['frames'] for v in records)==801
 report['cases'][label]=records
base=[sys.executable,str(Path(__file__).with_name('render_movies.py')),*map(str,cases),'--labels','Heated f-plane','Heated beta-plane','--title','Jupiter 55 N | identical initial eddies and 162 W/m² heating | days 0–100','--tracer','--fps','12','--crf','29']
for name,index in [('comparison',None),('f_plane',0),('beta_plane',1)]:
 target=out/f'{name}_0_100days.mp4';command=base+['--output',str(target)]
 if index is not None:command+=['--case-index',str(index)]
 subprocess.run(command,check=True)
 m=json.loads(target.with_suffix('.json').read_text());assert m['times_seconds']==[i*10800 for i in range(801)]
 # Publication copies retain all frames, dimensions and fixed physical scales.
 original=target.with_name(target.stem+'_original.mp4');shutil.copyfile(target,original)
 crf=31
 while target.stat().st_size>12_000_000 and crf<=39:
  temp=target.with_name(target.stem+'_web.mp4')
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
  temp.replace(target);m['web_reencode_crf']=crf;crf+=2
 assert target.stat().st_size<24_000_000
 count,duration=imageio_ffmpeg.count_frames_and_secs(str(target));assert count==801 and abs(duration-801/12)<.03
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(target),'-f','null','-'],check=True)
 m.update(final_mp4_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),full_decode_passed=True,compute_location='Unity')
 target.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
 report['media'][target.name]={'frames':count,'duration_seconds':duration,'bytes':target.stat().st_size,'sha256':m['final_mp4_sha256'],'full_decode_passed':True}
 print('VERIFIED',target.name,flush=True)
report['status']='COMPLETE';(out/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
print('READY_FOR_GITHUB_PUBLICATION',out,flush=True)

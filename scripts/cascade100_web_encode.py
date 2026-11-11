"""Encode smaller GitHub web copies on Unity; preserve dimensions and all 801 frames."""
from pathlib import Path
import shutil, json, subprocess, hashlib
import imageio_ffmpeg
src=Path('rendered/cascade100-full100-gallery');out=Path('rendered/cascade100-full100-web');out.mkdir(exist_ok=True)
for f in src.iterdir():
 if f.is_file():shutil.copyfile(f,out/f.name)
r=json.loads((out/'provenance.json').read_text())
for name,record in r['media'].items():
 f=out/name;m=json.loads(f.with_suffix('.json').read_text()); crf=35
 while f.stat().st_size>14000000:
  assert crf<=45
  temp=f.with_name(f.stem+'_tmp.mp4')
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-v','error','-i',str(src/name),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True)
  temp.replace(f);m['web_reencode_crf']=crf;crf+=2
 n,secs=imageio_ffmpeg.count_frames_and_secs(str(f));assert n==801 and abs(secs-801/12)<.03
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(f),'-f','null','-'],check=True)
 m['original_unity_mp4_sha256']=record['sha256'];m['final_mp4_sha256']=hashlib.sha256(f.read_bytes()).hexdigest();m['web_encoding_note']='Smaller H.264 web copy; unchanged dimensions, 801 frames, physical times and color scales.'
 f.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n');record.update(bytes=f.stat().st_size,sha256=m['final_mp4_sha256'],full_decode_passed=True);print('VERIFIED',name,record,flush=True)
r['web_encoder_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest();(out/'provenance.json').write_text(json.dumps(r,indent=2)+'\n')

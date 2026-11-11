"""Smaller web copies of the three long legacy movies; all saved frames retained."""
from pathlib import Path
import json,shutil,subprocess,hashlib
import imageio_ffmpeg
src=Path('rendered/legacy-campaigns-20261111');out=Path('rendered/legacy-campaigns-20261111-web');out.mkdir(exist_ok=True)
for f in src.iterdir():
 if f.is_file():shutil.copyfile(f,out/f.name)
for f in out.glob('*.mp4'):
 if f.stat().st_size<=14000000:continue
 m=json.loads(f.with_suffix('.json').read_text());original=m['sha256'];crf=28
 while f.stat().st_size>14000000:
  assert crf<=40
  temp=f.with_name(f.stem+'_tmp.mp4');subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-v','error','-i',str(src/f.name),'-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True);temp.replace(f);m['web_crf']=crf;crf+=2
 n,secs=imageio_ffmpeg.count_frames_and_secs(str(f));assert n==m['frames'];subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(f),'-f','null','-'],check=True)
 m.update(sha256=hashlib.sha256(f.read_bytes()).hexdigest(),original_unity_mp4_sha256=original,web_encoding_note='Unchanged dimensions, saved frames, physical times and color scale.');f.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n');print('VERIFIED',f.name,f.stat().st_size,flush=True)

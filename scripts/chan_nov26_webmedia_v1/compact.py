from pathlib import Path
import shutil,json,subprocess,hashlib,imageio_ffmpeg
src=Path('rendered/chan-nov26-anelastic10000');out=Path('rendered/chan-nov26-anelastic10000-web');out.mkdir(exist_ok=False)
media=json.loads((src/'media.json').read_text());assert media['full_decode_passed']
for p in src.iterdir():
 if p.suffix=='.mp4':
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-i',str(p),'-c:v','libx264','-preset','fast','-b:v','850k','-maxrate','1000k','-bufsize','2000k','-threads','4','-pix_fmt','yuv420p','-movflags','+faststart',str(out/p.name)],check=True)
  n,secs=imageio_ffmpeg.count_frames_and_secs(str(out/p.name));assert n==2001
  subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(out/p.name),'-f','null','-'],check=True)
  print(p.name,'COMPACT_FULL_DECODE_PASS',flush=True)
 elif p.name!='media.json':shutil.copy2(p,out/p.name)
media['web_encoding']={'bitrate_kbps':850,'resolution_unchanged':True,'frames_unchanged':True,'source_verified_render':str(src),'reason':'GitHub browser rejected larger movie upload; original renders retained on Unity.'}
media['assets']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'media.json').write_text(json.dumps(media,indent=2)+'\n');print('ALL_MEDIA_VERIFIED',flush=True)

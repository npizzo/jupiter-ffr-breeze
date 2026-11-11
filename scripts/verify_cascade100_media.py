from pathlib import Path
import subprocess,json,hashlib
import imageio_ffmpeg
p=Path('rendered/cascade100-full100-gallery');r=json.loads((p/'provenance.json').read_text());assert r['status']=='COMPLETE';checks=[]
for name,m in r['media'].items():
 f=p/name;meta=json.loads(f.with_suffix('.json').read_text());assert meta['frames']==801 and meta['times_seconds']==[i*10800 for i in range(801)]
 assert hashlib.sha256(f.read_bytes()).hexdigest()==m['sha256']
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(f),'-f','null','-'],check=True)
 checks.append({'file':name,'frames':801,'full_decode_passed':True,'start_day':0,'end_day':100})
(p/'media-validation.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','checks':checks},indent=2)+'\n');print('PASS',checks)

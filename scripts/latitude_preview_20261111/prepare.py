"""Freeze a stable saved-frame inventory and render an explicitly partial campaign on Unity."""
from pathlib import Path
import sys,json,tomllib,hashlib,subprocess,shutil,time,datetime
import numpy as np
import h5py,imageio_ffmpeg
root=Path.cwd();lat=int(sys.argv[1]);assert lat in (35,55)
source=root/('data/cascade35-v1/cascade35_heated' if lat==35 else 'data/cascadebeta100-v1/cascadebeta100_heated')
out=root/'rendered/latitude-preview-20261111'/f'latitude{lat}';out.mkdir(parents=True,exist_ok=False)
view=out/'view';view.mkdir()
frames=[];metadata=[];capture=time.time()
for segment in sorted(source.glob('segment_*'),key=lambda p:int(p.name.split('_')[-1])):
 mp=segment/'metadata.toml'
 if not mp.exists():continue
 raw=mp.read_bytes();m=tomllib.loads(raw.decode());dst=view/segment.name;dst.mkdir();(dst/'metadata.toml').write_bytes(raw)
 metadata.append({'segment':segment.name,'original_status':m['status'],'metadata_sha256':hashlib.sha256(raw).hexdigest(),'script_sha256':m['script_sha256'],'configuration_sha256':m['configuration_sha256']})
 shutil.copyfile(dst/'metadata.toml',out/f'{segment.name}_metadata.toml')
 for f in sorted(segment.glob('frame_*.jld2')):
  before=f.stat()
  if m['status']!='COMPLETE' and capture-before.st_mtime<30:continue
  with h5py.File(f) as q:
   t=float(q['t'][()]);v=q['vorticity_slices'][:];c=q['tracer_slices'][:];h=q['video_heights'][:]
  after=f.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'Frame changed during read'
  assert v.shape==c.shape==(3,m['ny'],m['nx']) and np.isfinite(v).all() and np.isfinite(c).all()
  assert np.array_equal(h,m['video_heights_actual_m'])
  (dst/f.name).symlink_to(f)
  digest=hashlib.sha256(v.tobytes()+c.tobytes()).hexdigest()
  frames.append({'file':str(f.relative_to(root)),'time_seconds':t,'bytes':after.st_size,'vorticity_tracer_native_bytes_sha256':digest})
assert frames and [f['time_seconds'] for f in frames]==[i*10800 for i in range(len(frames))]
last=frames[-1]['time_seconds']/86400
if lat==35:assert last>=20,'35 N publication must reach day 20'
manifest={'status':'FROZEN_PARTIAL_INVENTORY','latitude':lat,'captured_utc':datetime.datetime.fromtimestamp(capture,datetime.timezone.utc).isoformat(),'start_day':0,'end_day':last,'frames':frames,'segments':metadata,'compute_location':'Unity','original_fields_unchanged':True,'partial_campaign':True}
(out/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
movie=out/f'latitude{lat}_preview.mp4'
command=[sys.executable,str(Path(__file__).with_name('render_movies.py')),str(view),str(view),'--labels',f'{lat}°N heated β-plane',f'{lat}°N heated β-plane','--title',f'{lat}°N heated β-plane | days 0–{last:g} | partial campaign | fixed scales','--tracer','--case-index','0','--fps','12','--crf','29','--output',str(movie)]
subprocess.run(command,check=True)
m=json.loads(movie.with_suffix('.json').read_text());assert m['times_seconds']==[f['time_seconds'] for f in frames]
assert m['heights_m']==list(h) and m['frames']==len(frames)
original=movie.with_name(movie.stem+'_original.mp4');shutil.copyfile(movie,original)
for crf in [31,33,35,37]:
 if movie.stat().st_size<=12_000_000:break
 temp=movie.with_name(movie.stem+'_web.mp4')
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-i',str(original),'-an','-c:v','libx264','-crf',str(crf),'-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 temp.replace(movie);m['web_reencode_crf']=crf
assert movie.stat().st_size<24_000_000
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==len(frames) and abs(duration-count/12)<.04
subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
# Finished QA stills at start, midpoint, end; no rendering or decoding on the client.
qa=[]
for label,idx in [('start',0),('middle',count//2),('end',count-1)]:
 target=out/f'qa_{label}.jpg'
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-ss',str(idx/12),'-i',str(movie),'-frames:v','1','-q:v','2',str(target)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL);qa.append(target.name)
m.update(partial_campaign=True,start_day=0,end_day=last,compute_location='Unity',slurm_job_id=__import__('os').environ.get('SLURM_JOB_ID'),full_decode_passed=True,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),bytes=movie.stat().st_size,
         vorticity_scale_scope='Fixed independently for each depth over this movie; consult physical colorbar when comparing latitude movies.',
         tracer_center_basis='Same fixed center and width as the published full100 f-plane contrast movies, not fitted per frame or latitude.')
movie.with_suffix('.json').write_text(json.dumps(m,indent=2)+'\n')
names=[movie.name,movie.with_suffix('.png').name,movie.with_suffix('.json').name,'source_manifest.json']+[p.name for p in out.glob('segment_*_metadata.toml')]+qa
report={'status':'COMPLETE','latitude':lat,'end_day':last,'frames':count,'duration_seconds':duration,'movie':movie.name,'poster':movie.with_suffix('.png').name,'full_decode_passed':True,'compute_location':'Unity','assets':{n:{'bytes':(out/n).stat().st_size,'sha256':hashlib.sha256((out/n).read_bytes()).hexdigest(),'git_blob_sha1':hashlib.sha1(('blob '+str((out/n).stat().st_size)+'\0').encode()+(out/n).read_bytes()).hexdigest()} for n in names}}
(out/'PASS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='assets'}),flush=True)

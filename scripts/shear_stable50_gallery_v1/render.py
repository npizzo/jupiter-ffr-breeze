from pathlib import Path
import sys,subprocess,json,hashlib,tomllib
import imageio_ffmpeg
out=Path(sys.argv[1]);joined=out/'joined';media=out/'publication';media.mkdir(exist_ok=False)
cases=['stable_cyclonic','stable_none','stable_anticyclonic']
labels=['Cyclonic initial shear','No initial shear','Anticyclonic initial shear']
manifest={'status':'COMPLETE','simulation_days':[0,50],'frames':401,'fps':12,'movies':[],'source_simulations':{}}
for i,name in enumerate(cases):
 movie=media/(name+'_0_50.mp4')
 subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),
   *[str(joined/c) for c in cases],'--labels',*labels,'--title','Stable upper layer | dry Breeze f-plane, 55°N',
   '--tracer','--case-index',str(i),'--output',str(movie)],check=True)
 ffmpeg=imageio_ffmpeg.get_ffmpeg_exe()
 subprocess.run([ffmpeg,'-v','error','-threads','1','-i',str(movie),'-f','null','-'],check=True)
 # A decoded contact sheet stays on Unity; only finished publication media are transferred.
 subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-vf','fps=1/8,scale=600:-1,tile=2x2','-frames:v','1','-threads','1',str(media/(name+'_contact.jpg'))],check=True)
 metadata=tomllib.loads((joined/name/'metadata.toml').read_text())
 manifest['source_simulations'][name]={k:metadata[k] for k in ('status','simulated_seconds','nx','ny','nz','f','shear_rossby','heat_flux','cap_n2','video_heights_actual_m','script_sha256','configuration_sha256','gpu_name')}
 manifest['movies'].append(dict(case=name,file=movie.name,bytes=movie.stat().st_size,sha256=hashlib.sha256(movie.read_bytes()).hexdigest(),decode_check='PASS',duration_seconds=401/12))
 print('VERIFIED_MOVIE',name,movie.stat().st_size,flush=True)
manifest['assets']={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in media.iterdir() if p.is_file()}
(media/'publication.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('STABLE_THREE_MOVIES_COMPLETE',flush=True)

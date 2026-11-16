"""Unity publication: exact integrated work, explicit sampling limits, web MP4."""
from pathlib import Path
import json,sys,subprocess,hashlib,shutil
import numpy as np
import imageio_ffmpeg
old,new,render,out=map(Path,sys.argv[1:]);out.mkdir(exist_ok=False,parents=True)
r=json.loads((render/'report.json').read_text())
sign=np.genfromtxt(render/'vorticity_sign_history.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
for case,d in r['cases'].items():
    src=(old if case in ['forced','unforced'] else new)/case
    a=np.genfromtxt(src/'statistics.csv',delimiter=',',names=True)
    for name,w in d['windows'].items():
        lo,hi=map(float,name.split('_'))
        i=int(np.argmin(abs(a['growth_times']-lo)));j=int(np.argmin(abs(a['growth_times']-hi)))
        assert abs(a['growth_times'][i]-lo)<1e-6 and abs(a['growth_times'][j]-hi)<1e-6
        duration=a['t'][j]-a['t'][i]
        w['sampled_instantaneous_work']={k:w.pop(k) for k in ['mean_injection_work','mean_viscous_work','mean_drag_work']}
        w['integrated_work_rates']={k:float((a['energy_'+k][j]-a['energy_'+k][i])/duration) for k in ['advection','injection','viscosity','drag']}
        w['energy_storage_rate']=float((a['energy'][j]-a['energy'][i])/duration)
        w['integrated_balance_residual']=w['energy_storage_rate']-sum(w['integrated_work_rates'].values())
        assert abs(w['integrated_balance_residual'])<1e-10
    sel=(sign['case']==case)&(sign['growth_times']>=600-1e-6)
    d['late_mean_skewness']={k:float(sign[k][sel].mean()) for k in ['full_skewness','lowpass_skewness']}
    sel=(a['growth_times']>=600-1e-6)
    d['late_mean_low_k_energy_fraction']=float(a['low_k_energy_fraction'][sel].mean())
    d['flux_sampling_note']='Post-pulse samples every 1 time unit (5 pulse periods); not an average over pulse phase.' if case in ['random_sign','cyclonic_sign'] else 'Uniform 1-time-unit samples, with duplicate initial records removed.'
    d['window_rate_note']='Use integrated_work_rates for budget comparisons; sampled_instantaneous_work retains the original, phase-biased diagnostic for provenance.'
for name in ['poster.png','energy_transfer_history.png','spectra_flux.png','vorticity_sign_history.png','vorticity_sign_history.csv']:
    shutil.copyfile(render/name,out/name)
ff=imageio_ffmpeg.get_ffmpeg_exe();movie=out/'four_controls_120_1200.mp4'
subprocess.run([ff,'-v','error','-i',str(render/movie.name),'-vf','scale=1280:-2','-an','-c:v','libx264','-preset','medium','-threads','1','-b:v','2200k','-maxrate','2400k','-bufsize','4800k','-pix_fmt','yuv420p','-movflags','+faststart',str(movie)],check=True)
count,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));assert count==901 and movie.stat().st_size<25*1024**2
subprocess.run([ff,'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
r['status']='COMPLETE_REVIEWED_WITH_LIMITATIONS'
r['movie'].update(file=movie.name,frames=count,duration_seconds=duration,full_decode_passed=True,streaming_copy=True,width=1280)
r['original_movie_sha256']=hashlib.sha256((render/movie.name).read_bytes()).hexdigest()
r['publication_changes']=['Use exact time-integrated native energy work for window budgets; retain sampled values separately.', 'Stochastic flux and spectral samples remain explicitly post-pulse phase-conditioned.', 'All 901 video frames retained; browser copy scaled to 1280 pixels wide and encoded below GitHub browser-upload limit.']
r['assets']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'report.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
print(json.dumps({'movie_bytes':movie.stat().st_size,'frames':count,'cases':{c:{'late':d['windows']['600_1200'],'mean_skew':d['late_mean_skewness'],'mean_low_k_fraction':d['late_mean_low_k_energy_fraction']} for c,d in r['cases'].items()}},indent=2),flush=True)

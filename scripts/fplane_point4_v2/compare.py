"""Matched-case summary and synchronized movies; Unity only."""
from pathlib import Path
import sys,json,subprocess,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio_ffmpeg
root=Path(sys.argv[1]);results=Path(sys.argv[2]);out=results/'comparisons';out.mkdir(exist_ok=False)
names=(Path(__file__).with_name('case-order.txt')).read_text().splitlines()
reports={name:json.loads((results/name/'diagnostics.json').read_text()) for name in names}
for name in names:assert json.loads((results/name/'PASS.json').read_text())['full_decode_passed']
pairs=[('late_heated_dt15','late_heated_fine','Mature heated f-plane | 256² over 512² | days 90–100'),('early_replay_heated','early_replay_control','Original initial turbulence | heated over unheated | days 0–5'),('late_heated','late_control','Heated over unheated | days 90–100'),('rest_q162','rest_q162_fine','256² over 512² | from rest | 162 W/m²'),('rest_q162','rest_q162_flipf','Positive f over negative f | from rest | cyclones have ζ/f > 0'),('late_restart_control','late_restart_reversed','Unchanged over barotropic-reversed mature state | days 100–110')]
movies=[]
for a,b,title in pairs:
 movie=out/(a+'_over_'+b+'.mp4')
 subprocess.run([sys.executable,str(Path(__file__).with_name('render_movies.py')),str(root/a),str(root/b),'--labels',a,b,'--title',title,'--tracer','--fps','4','--crf','28','--output',str(movie)],check=True)
 frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie))
 subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True)
 assert movie.stat().st_size<24_000_000
 movies.append({'path':movie.name,'frames':frames,'duration':duration,'full_decode_passed':True})
summary={}
fig,ax=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
for name in names:
 r=reports[name];rows=r['snapshots'];t=np.array([x['time_days'] for x in rows]);s=[x['moments']['interior']['full_Z']['0']['raw_skewness'] for x in rows]
 tau=[x['source']['interior']['effective_tau_seconds_through_origin'] for x in rows]
 summary[name]={'grid':r['metadata']['nx'],'flux':r['metadata']['heat_flux'],'f':rows[0]['f'],'start_day':t[0],'end_day':t[-1],
  'final_skewness':s[-1],'final_core_fraction':rows[-1]['fits']['interior'].get('measured_core_fraction_zeta_over_rms_gt1'),
  'pooled_interior_relation':r['pooled_fits']['interior'],'tau_seconds':tau,
  'caution':'Short-window response; no sign-selection or third-invariant verdict from these numbers alone.'}
 if name.startswith('rest_q'):
  if name in ['rest_q162_flipf','rest_q162_fine']:continue
  ax[0,0].plot(t,s,label=f"{r['metadata']['heat_flux']:g} W/m²")
 if name in ['rest_q162','rest_q162_fine']:
  ax[0,1].plot(t,tau,label=name)
 if name in ['rest_q162','rest_q162_flipf']:
  ax[1,0].plot(t,s,label=name)
 if name in ['late_restart_control','late_restart_reversed']:
  ax[1,1].plot(t,s,label=name)
for a,title in zip(ax.ravel(),['Heat-flux sensitivity from rest','Effective τ: resolution sensitivity','Raw Z skewness after reversing f','Raw Z skewness after reversing mature flow']):
 a.set_title(title);a.set_xlabel('Earth days');a.legend(fontsize=8);a.axhline(0,color='.6',lw=.6)
ax[0,1].set_ylabel('τ [s]; fitted descriptive parameter')
for a in [ax[0,0],ax[1,0],ax[1,1]]:a.set_ylabel('S_Z (interior, density depth mean)')
fig.savefig(out/'campaign_comparison.png',dpi=150);plt.close(fig)
(out/'comparison.json').write_text(json.dumps({'compute_location':'Unity','cases':summary,'movies':movies},indent=2)+'\n')
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','assets':assets,'movies':movies},indent=2)+'\n')
print('POINT4_MATCHED_COMPARISONS_READY',out,flush=True)

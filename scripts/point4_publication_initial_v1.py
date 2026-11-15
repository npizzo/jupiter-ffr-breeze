"""Initial completed-case synthesis and media packaging. Run on Unity only."""
from pathlib import Path
import json,sys,shutil,subprocess,hashlib,datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio_ffmpeg
root=Path.cwd(); src=root/'rendered/fplane-point4-v2';out=root/'rendered/fplane-analysis-initial-v1'
out.mkdir(exist_ok=False)
names=[n for n in (root/'scripts/fplane_point4_v2/case-order.txt').read_text().splitlines() if n!='late_heated_fine']
reports={};series={};summary={}
for n in names:
 d=json.loads((src/n/'diagnostics.json').read_text());receipt=json.loads((src/n/'PASS.json').read_text());assert receipt['full_decode_passed']
 reports[n]=d;rows=[]
 for s in d['snapshots']:
  m=s['moments']['interior'];v=s['source']['interior'];fit=s['fits']['interior']
  rows.append({'day':s['time_days'],'f':s['f'],'Z_rms':m['full_Z']['0']['rms'],'raw_skewness':m['full_Z']['0']['raw_skewness'],'central_skewness':m['full_Z']['0']['central_skewness'],'nonzonal_raw_skewness':m['nonzonal_Z']['0']['raw_skewness'],'three_dimensional_raw_skewness':m['three_dimensional_uniform_volume']['raw_skewness'],'rotation_signed_3D_skewness':m['FuSun_reference']['rotation_signed_skewness_volume'],'reference_2p5_Ro':m['FuSun_reference']['2p5_Ro_volume'],'mean_Z':m['full_Z']['0']['mean'],'cyclonic_area_fraction':m['full_Z']['0']['cyclonic_fraction'],'tau_seconds':v['effective_tau_seconds_through_origin'],'tau_R2':v['tau_fit_R2'],'signed_V_active':v['rotation_signed_V_in_active_convection_mean'],'signed_H_ring':v['rotation_signed_H_in_200km_core_ring_mean'],'third_moment_sources':v['third_moment_source_s_inv4'],'quadratic_R2':fit.get('quadratic_R2'),'logistic_R2':fit.get('logistic_R2')})
 series[n]=rows
 bs=d['finite_interval_budgets'];rr=[b['unclosed_remainder_rms']/b['storage_rms'] for b in bs if b['storage_rms']>0]
 summary[n]={'grid':[d['metadata']['nx'],d['metadata']['ny'],d['metadata']['nz']],'heat_flux_W_m2':d['metadata']['heat_flux'],'first':rows[0],'last':rows[-1],'pooled_fits':d['pooled_fits'],'median_remainder_over_storage':float(np.median(rr)),'last_interval_budget':bs[-1],'interval_budgets':bs,'media_validation':receipt,'limitations':d['limitations']}
 for a,b in [(n+'.mp4',n+'.mp4'),(n+'.png',n+'.png'),('diagnostic_summary.png',n+'_diagnostics.png')]:shutil.copy2(src/n/a,out/b)
 # Public movie provenance; raw Unity data remain remote.
 provenance=json.loads((src/n/(n+'.json')).read_text());(out/(n+'_movie.json')).write_text(json.dumps(provenance,indent=2))
fig,ax=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
for n in ['rest_q5p4','rest_q54','rest_q162','rest_q540']:
 r=series[n];t=[x['day'] for x in r];label=f"{summary[n]['heat_flux_W_m2']:g} W/m²"
 ax[0,0].plot(t,[x['raw_skewness'] for x in r],label=label)
 ax[0,1].plot(t,[x['Z_rms'] for x in r],label=label)
for n in ['rest_q162','rest_q162_fine','rest_q162_flipf']:
 r=series[n];ax[1,0].plot([x['day'] for x in r],[None if x['raw_skewness'] is None else np.sign(x['f'])*x['raw_skewness'] for x in r],label=n.replace('rest_q',''))
for n in ['late_restart_control','late_restart_reversed']:
 r=series[n];label=n.replace('late_restart_','');ax[1,1].plot([x['day'] for x in r],[x['central_skewness'] for x in r],label=label+' central');ax[1,1].plot([x['day'] for x in r],[x['raw_skewness'] for x in r],ls='--',label=label+' raw')
for a,title,y in zip(ax.ravel(),['Heating sensitivity: from rest','Amplitude matters as well as skewness','Resolution and rotation sign at 162 W/m²','Reversal: raw and central moments disagree'],['Raw skewness of Z','RMS Z [s⁻¹]','sign(f) × raw skewness of Z','Skewness of Z']):
 a.set(title=title,xlabel='Earth days',ylabel=y);a.legend(fontsize=8);a.axhline(0,color='.6',lw=.5)
fig.suptitle('Initial f-plane analysis | interior, density-weighted depth-mean vorticity Z');fig.savefig(out/'initial_synthesis.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(1,3,figsize=(16,5),constrained_layout=True)
for n in ['late_heated','late_heated_dt15','late_control']:
 r=series[n];ax[0].plot([x['day'] for x in r],[x['central_skewness'] for x in r],label=n)
for n in ['late_heated','rest_q162','rest_q162_fine']:
 r=series[n];ax[1].plot([x['day'] for x in r],[x['signed_V_active'] for x in r],label=n+' V active');ax[1].plot([x['day'] for x in r],[x['signed_H_ring'] for x in r],ls='--',label=n+' H ring')
labels=['early_replay_heated','late_heated','rest_q162','rest_q162_fine','late_restart_reversed'];ax[2].bar(range(len(labels)),[summary[n]['median_remainder_over_storage'] for n in labels]);ax[2].set_xticks(range(len(labels)),['early heat','late heat','rest 256','rest 512','reversed'],rotation=30);ax[2].axhline(1,color='black',ls=':')
ax[0].set(title='Mature replay and timestep sensitivity',xlabel='Earth days',ylabel='Central skewness of Z');ax[0].legend(fontsize=8)
ax[1].set(title='Conditioned source signs (different masks)',xlabel='Earth days',ylabel='sign(f) × source mean [s⁻²]');ax[1].legend(fontsize=6)
ax[2].set(title='The centered budget is not closed',ylabel='Median residual RMS / storage RMS')
fig.savefig(out/'source_audit.png',dpi=150);plt.close(fig)
pairs=[('early_replay_heated','early_replay_control','Early turbulence | heated over zero bottom heat | days 0–5'),('late_heated','late_control','Mature turbulence | heated over zero bottom heat | days 90–100'),('rest_q162','rest_q162_fine','From rest at 162 W/m² | 256² over 512²'),('rest_q162','rest_q162_flipf','From rest | positive f over negative f'),('late_restart_control','late_restart_reversed','Mature state | unchanged over barotropic reversal | days 100–110')]
movies=[]
for a,b,title in pairs:
 n=a+'_over_'+b;movie=out/(n+'.mp4')
 subprocess.run([sys.executable,str(root/'scripts/fplane_point4_v2/render_movies.py'),str(root/'data/fplane-point4-v2'/a),str(root/'data/fplane-point4-v2'/b),'--labels',a,b,'--title',title,'--tracer','--fps','4','--crf','28','--output',str(movie)],check=True)
 frames,duration=imageio_ffmpeg.count_frames_and_secs(str(movie));subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-xerror','-i',str(movie),'-f','null','-'],check=True);assert movie.stat().st_size<24000000
 movies.append({'name':n,'title':title,'frames':frames,'duration_seconds':duration,'full_decode_passed':True})
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'compute_location':'Unity','completed_cases':names,'excluded_pending_case':'late_heated_fine','summary':summary,'series':series,'comparisons':movies,'interpretation_scope':'Initial 14-case analysis. No statistical ensemble, equilibrium, inverse-cascade flux, or native closed-budget claim.'}
# Remove absolute remote directory prefix from public metadata; no raw fields leave Unity.
(out/'analysis.json').write_text(json.dumps(report,indent=2).replace(str(root)+'/',''))
for p in out.glob('*.json'):p.write_text(p.read_text().replace(str(root)+'/',''))
assets={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.is_file()}
(out/'PASS.json').write_text(json.dumps({'status':'COMPLETE','compute_location':'Unity','movies':movies,'individual_movies_full_decode_passed':True,'assets':assets},indent=2))
print('INITIAL_ANALYSIS_PUBLICATION_READY',flush=True)

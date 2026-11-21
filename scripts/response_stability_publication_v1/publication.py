"""Read-only final campaign review and finished-asset bundle; Unity only."""
from pathlib import Path
import json,csv,shutil,hashlib,subprocess,tarfile
import numpy as np
import imageio_ffmpeg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('data/response-stability-publication-v1');root.mkdir(exist_ok=False)
def rows(p):
 with open(p) as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
response={'control':'2084628','plus':'2084629','minus':'2084630'}
stability={'control':'2084284','half':'2084407','full':'2084408'}
report={'response':{},'stability':{},'scope':'One parent, one location/template, finite-time experiments. No spontaneous formation or eigenmode claim.'}
for family,jobs in [('response',response),('stability',stability)]:
 review=Path('data/anelastic-response-v1-review' if family=='response' else 'data/formed-stability-v1-review')
 table=rows(review/'comparison/comparison.csv');summary=json.loads((review/'comparison/summary.json').read_text())
 report[family]['comparison']=summary;report[family]['initial_comparison']=table[0]
 shutil.copy(review/'comparison/comparison.csv',root/f'{family}-comparison.csv')
 shutil.copy(review/'comparison/response.png',root/f'{family}-response.png')
 for case,job in jobs.items():
  base=Path(f'data/anelastic-response-v1-{case}-{job}/parent' if family=='response' else f'data/formed-stability-v1-{case}-{job}/parent')
  # Completion receipts and actual final diagnostics agree.
  import tomllib
  done=tomllib.loads((base.parent/'COMPLETE.toml').read_text());assert done['status']=='COMPLETE'
  meta=tomllib.loads((base/'metadata.toml').read_text());d=rows(base/'diagnostics.csv')
  end=d[-1]['time'] if family=='response' else d[-1]['time_seconds']
  assert abs(end-(14000 if family=='response' else 102*86400))<1e-8
  vals={'job':job,'final':d[-1],'metadata':meta,'frame_metadata':json.loads((review/case/'movie.json').read_text())}
  if family=='response':
   rates=rows(base/'transfer.csv');vals['sample_count']=len(rates);vals['sampled_time_fraction']=sum(r['dt'] for r in rates)/2000
   keys=['energy_start','energy_end','energy_change','mean_flux_work','convective_stress_work','native_minus_centered_flux_work','coriolis_work','viscous_work','horizontal_forcing_work','closure_residual']
   vals['transfer_windows']=[]
   for lo,hi in [(12050,12500),(12500,13000),(13000,13500),(13500,14000),(12050,14000)]:
    selected=[r for r in rates if lo<=r['time']<=hi];assert selected
    means={k:float(np.mean([r[k]/r['dt'] for r in selected])) for k in keys if k not in ['energy_start','energy_end']}
    conv=np.array([r['convective_stress_work']/r['dt'] for r in selected]);unc=np.array([r['native_minus_centered_flux_work']/r['dt'] for r in selected]);err=np.array([r['closure_residual']/r['dt'] for r in selected]);scale=max(float(np.mean(abs(conv))),1e-25)
    vals['transfer_windows'].append({'start':lo,'end':hi,'samples':len(selected),'mean_sampled_rates':means,'mean_abs_convective_rate':scale,'mean_abs_native_discrepancy_over_convective':float(np.mean(abs(unc))/scale),'max_abs_closure_rate_over_convective':float(max(abs(err))/scale)})
   vals['initial_seed']=tomllib.loads((base/'PICKUP_PASS.toml').read_text())
   shutil.copy(base/'transfer.csv',root/f'{family}-{case}-transfer.csv')
  else:vals['initial_perturbation']=tomllib.loads((base/'PERTURBATION_PASS.toml').read_text())
  report[family][case]=vals
  shutil.copy(base/'diagnostics.csv',root/f'{family}-{case}-diagnostics.csv')
  movie=root/f'{family}-{case}.mp4';shutil.copy(review/case/'evolution.mp4',movie)
  shutil.copy(review/case/'movie.json',root/f'{family}-{case}-movie.json')
  ffmpeg=imageio_ffmpeg.get_ffmpeg_exe();duration=vals['frame_metadata']['frame_count']/10
  # Poster sampled from the validated movie; the last half second is available in every clip.
  subprocess.run([ffmpeg,'-v','error','-ss',str(max(duration-.5,0)),'-i',str(movie),'-frames:v','1','-q:v','2',str(root/f'{family}-{case}.jpg')],check=True)
  assert movie.stat().st_size<24_000_000
 if family=='stability':
  init=table[0];end=table[-1]
  report[family]['difference_energy_ratio_final_initial']={str(j):end[f'large_scale_difference_energy_{j}']/init[f'large_scale_difference_energy_{j}'] for j in [1,2]}
  report[family]['difference_amplitude_ratio_final_initial']={str(j):float(np.sqrt(end[f'large_scale_difference_energy_{j}']/init[f'large_scale_difference_energy_{j}'])) for j in [1,2]}
  report[family]['max_linearity_error']=max(r['half_full_linearity_error'] for r in table)
 else:
  report[family]['mode_windows']=[]
  for lo,hi in [(12005,12050),(12050,12500),(12500,13000),(13000,13500),(13500,14000)]:
   sel=[r for r in table if lo<=r['time']<=hi]
   report[family]['mode_windows'].append({'start':lo,'end':hi,'means':{k:float(np.mean([r[k] for r in sel])) for k in sel[0] if k!='time'}})
# Plot actual sampled work rates; these are not an integrated full-interval energy budget.
fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
for ax,case in zip(axs,response):
 rr=rows(root/f'response-{case}-transfer.csv')
 for key,label in [('convective_stress_work','Resolved convection'),('mean_flux_work','Mean flow'),('native_minus_centered_flux_work','Native − centered flux'),('viscous_work','Viscosity')]:
  ax.plot([r['time'] for r in rr],[r[key]/r['dt'] for r in rr],label=label,lw=.9)
 ax.set_title(case);ax.set_xlabel('Nondimensional time');ax.grid(alpha=.2)
axs[0].set_ylabel('Sampled low-mode energy work rate');axs[-1].legend(fontsize=8)
fig.savefig(root/'sampled-transfer.png',dpi=140);plt.close(fig)
(root/'review.json').write_text(json.dumps(report,indent=2)+'\n')
manifest={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in root.iterdir() if p.is_file()}
(root/'assets.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open(str(root)+'.tar','w') as tar:
 for p in root.iterdir():tar.add(p,arcname=p.name)
brief={'response_mode_windows':report['response']['mode_windows'],'response_transfer':{c:report['response'][c]['transfer_windows'][-1] for c in response},'stability_difference_amplitude_ratio':report['stability']['difference_amplitude_ratio_final_initial'],'stability_max_linearity_error':report['stability']['max_linearity_error'],'movie_bytes':{k:v['bytes'] for k,v in manifest.items() if k.endswith('.mp4')}}
print(json.dumps(brief,indent=2),flush=True)

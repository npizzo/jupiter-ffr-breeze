"""Unity-only budget comparison at native and common coarse resolution."""
from pathlib import Path
import json,sys,tomllib,hashlib
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
paths={'coarse_dt30':Path('data/ferrel35-face-v3/production'),'coarse_dt15':root/'coarse','fine_dt15':root/'fine'}
windows={'full_24h':(300,301),'adjustment_6h':(300,300.25),'later_18h':(300.25,301)}
report={'status':'PASS','compute_location':'Unity','windows_earth_days':windows,'cases':{},'interpretation':['Horizontal refinement only; vertical grid, physical coefficients and boundaries remain fixed.','The 15-second coarse run separates timestep-limit sensitivity from horizontal resolution sensitivity.','The initial fine state is a conservative divergence-preserving prolongation of the coarse day-300 state; it is not an independently equilibrated fine simulation.','Native-grid and common-grid norms are both reported; common-grid fine tendencies are averaged over pairs of meridional cells before measuring.','A smaller reconstruction/eddy ratio alone does not establish convergence; compare absolute contributions and early versus later windows.','WENO-minus-linear-product is a specified reconstruction correction, not automatically dissipative stress.','The column enthalpy diagnostic uses reference Exner cell-center weights, not a closed total-energy budget.']}
curves={}
for label,path in paths.items():
 meta=tomllib.loads((path/'metadata.toml').read_text());assert meta['status']=='COMPLETE' and meta['simulated_seconds']==301*86400
 files={}
 for p in sorted(path.glob('frame_*.jld2')):
  with h5py.File(p) as d:files[float(d['t'][()])/86400]=p
 assert all(t in files for t in [300.25,301])
 with h5py.File(files[301]) as d:
  rho=d['density'][:];y=d['y'][:];z=d['z'][:];pressure=d['pressure'][:]
  origin_heat=d['face_heat_origin'][:]
  assert float(d['face_origin_time'][()])==300*86400
  max_heat_closure=float(abs(d['face_heat_residual'][:]).max()/86400)
  assert max_heat_closure<1e-12
  steps=int(d['face_steps'][()]);assert steps==meta['native_budget_checks']
  full_face=d['face_integrated_momentum'][:];full_native=d['native_integrated_parts'][:]
  exact_error=float(np.max(abs(full_face.sum(0)-full_native[0]))/86400)
  assert exact_error<1e-12
 exner=(pressure/1e5)**((8.314462618/meta['molar_mass'])/meta['heat_capacity'])
 energy=meta['heat_capacity']*exner[:,None];dz=meta['height']/meta['nz']
 column=lambda q:float(dz*np.mean(np.sum(energy*q,axis=0)))
 def state(t):
  if t==300:
   return np.zeros_like(full_face),np.zeros_like(full_native),np.zeros((5,len(z),len(y))),origin_heat
  with h5py.File(files[t]) as d:
   return d['face_integrated_momentum'][:],d['native_integrated_parts'][:],d['face_integrated_heat'][:],d['face_heat_current'][:]
 def metrics(mom,native,yy):
  mask=(abs(yy)[None,:]<=2400000)&(z[:,None]>=10000)&(z[:,None]<=75000)
  weights=np.broadcast_to(rho[:,None],mask.shape)[mask].copy();weights/=weights.sum()
  rms=lambda q:float(np.sqrt(np.sum(weights*q[mask]**2)))
  mean=mom[0]+mom[1];eddy=mom[2]+mom[3];correction=mom[4]+mom[5]
  return {'mean_rms':rms(mean),'eddy_rms':rms(eddy),'reconstruction_rms':rms(correction),'reconstruction_to_eddy':rms(correction)/rms(eddy),
    'eddy_y_rms':rms(mom[2]),'eddy_z_rms':rms(mom[3]),'reconstruction_y_rms':rms(mom[4]),'reconstruction_z_rms':rms(mom[5]),
    'coriolis_rms':rms(native[1]),'diffusion_rms':rms(native[2]),'native_advection_rms':rms(native[0]),
    'face_sum_vs_native_rms':rms(mom.sum(0)-native[0])}
 rec={'grid':[meta['nx'],meta['ny'],meta['nz']],'dx_km':meta['lx']/meta['nx']/1000,'max_dt':meta['max_dt'],'native_steps':steps,'maximum_heat_closure_rho_theta_per_second':max_heat_closure,'maximum_momentum_partition_error_rho_u_per_second':exact_error,'source':str(path),'metadata_sha256':hashlib.sha256((path/'metadata.toml').read_bytes()).hexdigest(),'windows':{}}
 for window,(a,b) in windows.items():
  a0,b0=state(a),state(b);duration=(b-a)*86400
  mom=(b0[0]-a0[0])/duration/rho[None,:,None];native=(b0[1]-a0[1])/duration/rho[None,:,None]
  heat=(b0[2]-a0[2])/duration;storage=(b0[3]-a0[3])/duration
  mm,nn,yy=mom,native,y
  if label=='fine_dt15':
   mm=mom.reshape(6,len(z),256,2).mean(-1);nn=native.reshape(4,len(z),256,2).mean(-1);yy=y.reshape(256,2).mean(-1)
  rec['windows'][window]={'native_grid':metrics(mom,native,y),'common_coarse_grid':metrics(mm,nn,yy),'column_enthalpy_W_m2':dict(zip(['advection','diffusion','stable_cap','thermal_feedback','bottom_flux'],map(column,heat)))}
  rec['windows'][window]['column_enthalpy_W_m2']['storage']=column(storage)
  rec['windows'][window]['column_enthalpy_W_m2']['residual']=column(storage-heat.sum(0))
  if window=='later_18h':curves[label]=(yy,z,mm[2]+mm[3],mm[4]+mm[5])
 if (path/'remap_validation.toml').exists():rec['remap_validation']=tomllib.loads((path/'remap_validation.toml').read_text())
 report['cases'][label]=rec
report['preflight']=tomllib.loads((root/'preflight/PASS.toml').read_text())
# Keep source hashes, bounds, and actual remap energy changes in the public record.
report['analysis_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(out/'resolution_review.json').write_text(json.dumps(report,indent=2)+'\n')
fig,axs=plt.subplots(1,3,figsize=(15,4.5),constrained_layout=True)
labels=['256², 30 s','256², 15 s','512², 15 s'];colors=['#8996a7','#1f6f8b','#cc6b27']
for ax,metric,title in zip(axs,['reconstruction_to_eddy','eddy_rms','reconstruction_rms'],['Reconstruction / eddy RMS','Resolved eddy acceleration RMS (m/s²)','Reconstruction acceleration RMS (m/s²)']):
 for i,case in enumerate(paths):
  vals=[report['cases'][case]['windows'][w]['common_coarse_grid'][metric] for w in windows]
  ax.bar(np.arange(3)+(i-1)*.24,vals,.24,label=labels[i],color=colors[i])
 ax.set_xticks(range(3),['Full day','First 6 h','Final 18 h']);ax.set_title(title);ax.grid(axis='y',alpha=.15)
axs[0].legend(fontsize=8);fig.suptitle('35°N day-300 restart: horizontal resolution and timestep sensitivity\nAll tendencies compared on common coarse meridional cells; one day is not a convergence proof')
fig.savefig(out/'resolution_summary.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
for row,index in enumerate([2,3]):
 lim=max(np.max(abs(curves[k][index])[(abs(curves[k][0])[None,:]<=2400000)&(curves[k][1][:,None]>=10000)&(curves[k][1][:,None]<=75000)]) for k in paths)
 for col,k in enumerate(paths):
  yy,zz,_,_=curves[k];im=axs[row,col].pcolormesh(yy/1e6,zz/1000,curves[k][index],cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
  axs[row,col].set(title=labels[col],xlabel='Northward distance (1000 km)',ylabel='Height (km)',xlim=(-2.4,2.4),ylim=(10,75));fig.colorbar(im,ax=axs[row,col],label='m/s²')
fig.suptitle('Final 18 hours: resolved eddy acceleration (top), reconstruction correction (bottom)\nCommon coarse-grid averages; fixed shared ranges within each row')
fig.savefig(out/'resolution_budgets.png',dpi=130);plt.close(fig)
print(json.dumps({k:v['windows'] for k,v in report['cases'].items()},indent=2),flush=True)

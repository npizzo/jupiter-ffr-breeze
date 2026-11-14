"""Summarize the diagnostic branch on Unity; no simulation mutation."""
from pathlib import Path
import json,sys,tomllib,hashlib
import numpy as np,h5py
root=Path(sys.argv[1]);meta=tomllib.loads((root/'metadata.toml').read_text());files=sorted(root.glob('frame_*.jld2'))
assert meta['status']=='COMPLETE'
with h5py.File(files[-1]) as d:
 t0=float(d['native_budget_origin_time'][()]);t1=float(d['native_budget_end_time'][()]);duration=t1-t0
 assert duration>0 and int(d['native_budget_steps'][()])==meta['native_budget_checks']
 rho=d['density'][:];y=d['y'][:];z=d['z'][:]
 # HDF5 reverses Julia dimensions: (term,z,y).
 parts=d['native_integrated_parts'][:]/duration/rho[None,:,None]
 residual=d['native_integrated_residual'][:]/duration/rho[:,None]
 assert parts.shape==(4,len(z),len(y))
 mask=(abs(y[None,:])<=.4*meta['ly'])&(z[:,None]>=10000)&(z[:,None]<=75000)
 weights=np.broadcast_to(rho[:,None],mask.shape)[mask];weights/=weights.sum()
 rms=lambda a:float(np.sqrt(np.sum(weights*a[mask]**2)))
 values={k:rms(a) for k,a in zip(['advection','coriolis','diffusion','external_forcing'],parts)}
 values['closed_budget_residual']=rms(residual)
 assert np.isfinite(parts).all() and np.isfinite(residual).all()
 assert float(d['native_integrated_relative_residual'][()])<1e-5
names=['centered_mean_advection_y','centered_mean_advection_z','centered_eddy_advection_y','centered_eddy_advection_z','centered_native_advection_remainder']
times=[];records={k:[] for k in names};native_snap=[]
for p in files:
 with h5py.File(p) as d:
  times.append(float(d['t'][()]));native_snap.append(d['native_rhou_rhs'][:])
  for k in names:records[k].append(d[k][:])
times=np.array(times);assert times[0]>=t0 and times[-1]==t1 and np.all(np.diff(times)>0)
# Restart callbacks omit the duplicate checkpoint-time frame. Compare snapshots
# with native integrals over their actual common interval, not the full branch.
sample_duration=times[-1]-times[0]
with h5py.File(files[0]) as d:
 first_integrated=d['native_integrated_parts'][:]
with h5py.File(files[-1]) as d:
 native_common=(d['native_integrated_parts'][:]-first_integrated)/sample_duration/rho[None,:,None]
means={k:np.trapezoid(np.asarray(a),x=times,axis=0)/sample_duration/rho[:,None] for k,a in records.items()}
ns=np.trapezoid(np.asarray(native_snap),x=times,axis=0)/sample_duration/rho[None,:,None]
report={'status':'PASS','scope':'Numerical gate only; a closed total budget does not establish eddy-driven Ferrel cells','compute_location':'Unity','start_day':t0/86400,'end_day':t1/86400,'snapshot_mean_start_day':times[0]/86400,'snapshot_mean_end_day':times[-1]/86400,'native_steps':meta['native_budget_checks'],'native_time_integrated_acceleration_rms_m_s2':values,'centered_3hour_sampled_acceleration_rms_m_s2':{k:rms(v) for k,v in means.items()},'three_hour_sampling_difference_rms_m_s2':{k:rms(ns[i]-native_common[i]) for i,k in enumerate(['advection','coriolis','diffusion','external_forcing'])},'limitations':['Every-step native budget uses exact SSPRK3 weights; centered mean/eddy estimates use trapezoidal three-hour snapshots.','Native-minus-centered advection difference contains interpolation and numerical-transport effects; it is not automatically a physical eddy stress.','This gate validates diagnostic operation and conservation, not cell existence, resolution convergence or statistical equilibrium.'],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(root/'native_budget_review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

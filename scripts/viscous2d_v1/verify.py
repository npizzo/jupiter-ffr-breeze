"""Unity-only gate: actual unforced upscale transfer plus controlled forcing probes."""
from pathlib import Path
import sys,json,tomllib,hashlib
import numpy as np
root=Path(sys.argv[1]);results={};events={};initial_hashes=[]
for case in ['unforced','forced','random_sign','cyclonic_sign']:
 p=root/case;c=tomllib.loads((p/'config.toml').read_text());r=tomllib.loads((p/'review.toml').read_text())
 a=np.genfromtxt(p/'statistics.csv',delimiter=',',names=True)
 assert r['status']=='COMPLETED' and r['max_cfl']<.5
 assert c['cubic_velocity_drag']==0 and c['viscosity']==.0001
 assert np.all(np.isfinite(a.view(float)))
 assert abs(r['energy_increments'][3])<1e-20
 assert abs(a['energy'][-1]-a['energy'][0]-sum(r['energy_increments']))/a['energy'][0]<1e-8
 initial_hashes.append(hashlib.sha256((p/'zeta_0000000.f32').read_bytes()).hexdigest())
 S=np.genfromtxt(p/'spectra.csv',delimiter=',',names=True);nk=len(np.unique(S['k']));D={n:S[n].reshape(-1,nk) for n in S.dtype.names};k=D['k'][0]
 it=D['iteration'][:,0].astype(int);assert np.array_equal(it[it>0]%80,(81*np.arange(1,len(it[it>0])+1))%80)
 results[case]={'review':r,'config':c,'sample_count':len(it)}
 if case=='unforced':
  late=D['t'][:,0]>=4;flux=D['Pi_advection'][late]
  meanflux=[float(flux[:,k==K].mean()) for K in [4,6,8]]
  negative=[float((flux[:,k==K]<0).mean()) for K in [4,6,8]]
  elo=D['energy'][:,k<=4].sum(axis=1)
  transfer={'late_mean_flux_K4_K6_K8':meanflux,'negative_fraction_K4_K6_K8':negative,'low_k_energy_initial_final':[float(elo[0]),float(elo[-1])],'low_k_gain_over_initial_energy':float((elo[-1]-elo[0])/a['energy'][0])}
  assert max(meanflux)<0,transfer
  assert min(negative)>=.75,transfer
  assert transfer['low_k_gain_over_initial_energy']>=.01,transfer
  n=c['grid'][0];kx=np.fft.fftfreq(n)*n;k2=kx[:,None]**2+kx[None,:]**2;high=k2>=(n/3)**2
  tails=[]
  for f in sorted(p.glob('zeta_*.f32')):
   q=np.fromfile(f,dtype=np.float32).reshape((n,n),order='F').astype(float)
   power=abs(np.fft.fft2(q))**2;tails.append(float(power[high].sum()/power.sum()))
  transfer['max_enstrophy_fraction_at_k_ge_N_over_3']=max(tails)
  assert max(tails)<.01,transfer
  assert abs(r['energy_increments'][1])<1e-20
  results[case]['cascade_gate']=transfer
 elif case in ['random_sign','cyclonic_sign']:
  e=np.genfromtxt(p/'events.csv',delimiter=',',names=True);events[case]=e
  assert len(e)==10
  assert np.allclose(e['energy_quadratic'],c['expected_energy_injection_rate']*.2,rtol=1e-12)
  assert abs(e['energy_total'].sum()-r['energy_increments'][1])<1e-10
  assert np.all(e['core_sign']==1) if case=='cyclonic_sign' else set(e['core_sign'])=={-1,1}
 else:assert np.isclose(c['amplitude'],.0002)
assert len(set(initial_hashes))==1
for k in ['shift_x','shift_y','energy_quadratic','iteration']:assert np.array_equal(events['random_sign'][k],events['cyclonic_sign'][k])
report={'status':'PASS','compute_location':'Unity','scope':'Transient decaying inverse transfer; no stationary inertial-range or convergence claim','cases':results}
(root/'PASS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':'PASS','cascade_gate':results['unforced']['cascade_gate']},indent=2),flush=True)

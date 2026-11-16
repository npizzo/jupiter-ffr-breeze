"""Unity-only gate on the actual production driver, including pulse-phase coverage."""
from pathlib import Path
import json,sys,tomllib,hashlib
import numpy as np
root=Path(sys.argv[1]);configs={};events={};records={}
for case in ['forced','random_sign','cyclonic_sign']:
 p=root/case;c=tomllib.loads((p/'config.toml').read_text());r=tomllib.loads((p/'review.toml').read_text())
 a=np.genfromtxt(p/'statistics.csv',delimiter=',',names=True);s=np.genfromtxt(p/'spectra.csv',delimiter=',',names=True)
 assert r['status']=='COMPLETED' and r['max_cfl']<.5
 assert c['forcing_scale']==.2 and c['viscosity']==.001 and c['cubic_velocity_drag']==1
 assert np.all(np.isfinite(a.view(float)))
 assert abs(a['energy'][-1]-a['energy'][0]-sum(r['energy_increments']))/a['energy'][0]<1e-8
 its=np.unique(s['iteration']).astype(int);its=its[its>0]
 phases=its%40
 assert len(its)==40 and set(phases)==set(range(40)),(case,len(its),phases)
 assert np.array_equal(a['iteration'][1:].astype(int),its)
 configs[case]=c
 records[case]={'review':r,'sample_count':len(its),'distinct_phases':len(set(phases)),
 'initial_field_sha256':hashlib.sha256((p/'zeta_0000000.f32').read_bytes()).hexdigest()}
 if case!='forced':
  e=np.genfromtxt(p/'events.csv',delimiter=',',names=True);events[case]=e
  assert len(e)==201
  assert np.allclose(e['energy_quadratic'],c['expected_energy_injection_rate']*.2,rtol=1e-12)
  assert np.all(e['core_sign']==1) if case=='cyclonic_sign' else set(e['core_sign'])=={-1,1}
  assert abs(e['energy_total'].sum()-r['energy_increments'][1])<1e-10
 else:assert r['pulse_events']==0 and c['amplitude']==.001875
assert len({c['source_checkpoint_sha256'] for c in configs.values()})==1
assert len({v['initial_field_sha256'] for v in records.values()})==1
for k in ['shift_x','shift_y','energy_quadratic','iteration']:
 assert np.array_equal(events['random_sign'][k],events['cyclonic_sign'][k])
(root/'PASS.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','cases':records},indent=2)+'\n')
print('PASS: actual driver, three matched probes, all 40 pulse phases, integrated budgets and paired random streams')

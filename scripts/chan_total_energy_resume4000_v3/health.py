"""Endpoint release gate; no simulation modifications or external dependencies."""
from pathlib import Path
import subprocess,csv,hashlib,json,math,os,shutil,sys,tomllib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def check(receipt):
 c=tomllib.loads((HERE/'config.toml').read_text()); donor=ROOT/c['donor_root']; parent=donor/'parent'
 state=subprocess.check_output(['sacct','-j','2261539','-X','--noheader','--parsable2','--format=State'],text=True).strip()
 assert state.startswith('CANCELLED'), ('DONOR_NOT_STOPPED',state)
 assert not (parent/'GUARD_STOP.toml').exists(), 'DONOR_GUARD_STOP'
 meta=tomllib.loads((parent/'metadata.toml').read_text())
 # Donor was explicitly cancelled by the user; resume its last complete checkpoint.
 for k in ['nx','ny','nz','omega','beta','jet_speed','bottom_flux','cooling_time','heat','seed','perturbation','initial_rms','closure','advection_order']:
  assert meta[k]==c[k], ('PHYSICS_MISMATCH',k)
 rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((parent/'diagnostics.csv').open())]
 rows=[r for r in rows if r['time']<=c['restart_time']+1e-8]
 assert rows and abs(rows[-1]['time']-c['restart_time'])<1e-8
 assert all(all(math.isfinite(v) for v in r.values()) for r in rows)
 assert all(b['time']>a['time'] for a,b in zip(rows,rows[1:]))
 for r in rows:
  assert r['min_density']>0 and r['min_temperature']>0 and r['dt']>0
  assert r['mass_error']<c['mass_tolerance'] and r['constraint']<c['constraint_tolerance']
  assert r['max_mach']<c['anelastic_mach_stop' if c['equations']=='anelastic' else 'compressible_mach_stop']
 if c.get('formulation')=='TotalEnergy':
  assert tomllib.loads((parent/'FULL_GRID_ENERGY_PASS.toml').read_text())['status']=='PASS'
  assert all(abs(r['defect_over_bottom_heat'])<c['budget_relative_tolerance'] for r in rows if r['time']-rows[0]['time']>=c['budget_gate_start'])
  assert all(abs(r['reconstruction_minus_native'])/max(abs(r['energy_proxy']),1)<1e-10 for r in rows)
 cps=list((parent/'checkpoints').glob('checkpoint_iteration*.jld2'));assert len(cps)==1, 'AMBIGUOUS_CHECKPOINT'
 cp=cps[0];assert int(cp.stem.split('iteration')[-1])==int(rows[-1]['iteration']), 'CHECKPOINT_ITERATION_MISMATCH'
 assert cp.stat().st_size>1000000
 assert shutil.disk_usage(ROOT).free>15*1024**3, 'LESS_THAN_15_GIB_FREE'
 out=Path(receipt);out.parent.mkdir(parents=True,exist_ok=True)
 result={'status':'PASS','restart_time':c['restart_time'],'restart_checkpoint':str(cp.relative_to(ROOT)),
 'checkpoint_sha256':digest(cp),'diagnostics_sha256':digest(parent/'diagnostics.csv'),
 'donor_metadata_sha256':digest(parent/'metadata.toml'),'target_time':c['stop_time'],
 'energy_proxy_is_not_anelastic_invariant':c['equations']=='anelastic',
 'flagged_original_theta_diagnostic':c.get('flagged_diagnostic',False),
 'endpoint':rows[-1], 'peak_mach':max(r['max_mach'] for r in rows),
 'cap_and_profile_are_reported_not_equilibrium_gates':True}
 out.write_text(json.dumps(result,indent=2)+'\n')
 # Native driver reads this simple TOML; the JSON records the complete assessment.
 vals={k:result[k] for k in ['status','restart_time','restart_checkpoint','checkpoint_sha256','target_time']}
 vals.update({'donor_'+k:rows[-1][k] for k in ['kinetic_energy','energy_proxy','min_density','min_temperature']})
 out.with_suffix('.toml').write_text('\n'.join(k+' = '+json.dumps(v) for k,v in vals.items())+'\n')
 print(json.dumps(result),flush=True)
def verify(receipt):
 p=Path(receipt);r=json.loads(p.read_text());c=tomllib.loads((HERE/'config.toml').read_text())
 assert r['status']=='PASS' and r['restart_time']==c['restart_time'] and r['target_time']==c['stop_time']
 assert digest(ROOT/r['restart_checkpoint'])==r['checkpoint_sha256']
 assert digest(ROOT/c['donor_root']/'parent/diagnostics.csv')==r['diagnostics_sha256']
 assert digest(ROOT/c['donor_root']/'parent/metadata.toml')==r['donor_metadata_sha256']
 t=tomllib.loads(p.with_suffix('.toml').read_text())
 for k in ['status','restart_time','restart_checkpoint','checkpoint_sha256','target_time']:assert t[k]==r[k]
 for k in ['kinetic_energy','energy_proxy','min_density','min_temperature']:assert t['donor_'+k]==r['endpoint'][k]
 print('DONOR_HEALTH_AND_CHECKPOINT_VERIFIED',flush=True)
if __name__=='__main__':
 try:
  {'check':check,'verify':verify}[sys.argv[1]](sys.argv[2])
 except Exception as e:
  Path(sys.argv[2]).with_suffix('.failure.json').write_text(json.dumps({'status':'BLOCKED','reason':repr(e)},indent=2)+'\n')
  raise

# Run only on Unity; reads scalar diagnostics and hashes checkpoints in place.
import csv,json,sys,hashlib,tomllib,statistics
from pathlib import Path
root=Path(sys.argv[1]);receipt_path=Path(sys.argv[2]);receipt=json.loads(receipt_path.read_text());assert receipt["status"]=="PASS";report={};checks=[];lines=['status = "PASS"']
for label in ['case3','case4']:
 d=root/label
 meta=tomllib.loads((d/'metadata.toml').read_text())
 init=tomllib.loads((d/'INITIALIZATION.toml').read_text())
 assert meta['model_sha256']==receipt['manifest']['files']['scripts/chan_eddies_v3/model.jl']
 assert (meta['nx'],meta['ny'],meta['nz'])==(256,256,80)
 assert meta['omega']==({'case3':1.5,'case4':3.}[label]) and meta['initial_rms']==.136
 rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader((d/'diagnostics.csv').open())]
 r0=rows[0];last=rows[-1]
 maxres=max(abs(r['energy_residual_over_input']) for r in rows)
 masserr=max(abs(r['mass']/r0['mass']-1) for r in rows)
 ratios=[r['barotropic_ke']/r0['barotropic_ke'] for r in rows]
 timings=[(b['wall_seconds']-a['wall_seconds'])/(b['time']-a['time']) for a,b in zip(rows,rows[1:]) if a['time']>=2 and b['time']>a['time']]
 rate=statistics.median(timings)
 integral=sum(.5*(a['external_power']+b['external_power'])*(b['time']-a['time']) for a,b in zip(rows,rows[1:]))
 integrated=(last['total_energy']-r0['total_energy']-integral)/(16*meta['bottom_flux']*(last['time']-r0['time']))
 cp=list((d/'checkpoints').glob('*.jld2'));assert len(cp)==1
 h=hashlib.sha256()
 with cp[0].open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 target=meta['initial_turnovers']*init['turnover_time']
 ok=init['N2_cap_min']>0 and max(abs(init['N2_deep_min']),abs(init['N2_deep_max']))<.2 and meta['status']=='COMPLETE' and last['time']>=5-1e-8 and maxres<.1 and masserr<1e-7 and max(r['max_mach'] for r in rows)<.5 and min(ratios)>.5 and max(ratios)<2 and abs(integrated)<.1
 checks.append(ok)
 report[label]=dict(passed=ok,final_time=last['time'],max_energy_residual_over_input=maxres,integrated_energy_residual_over_input=integrated,max_mass_error=masserr,barotropic_ke_ratio_min=min(ratios),barotropic_ke_ratio_max=max(ratios),median_seconds_per_time_unit=rate,conditional_remaining_seconds=rate*(target-5),allocation_wall_limit_seconds=5400,integration_wall_limit_seconds=5100,full_target_within_first_allocation=rate*(target-5)<4800,initialization=init)
 lines += [label+'_checkpoint = '+json.dumps(str(cp[0].resolve())),label+'_checkpoint_sha256 = '+json.dumps(h.hexdigest())]
for key in ['mass','kinetic_energy']:
 checks.append(abs(report['case3']['initialization'][key]/report['case4']['initialization'][key]-1)<1e-10)
report['status']='PASS' if all(checks) else 'FAIL';report['scope']='Short numerical/adjustment screen only; no thermal equilibrium or long-time validity claim.'
(root/'SCREEN_REVIEW.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if not all(checks):sys.exit(2)
lines += ['cpu_receipt_sha256 = '+json.dumps(hashlib.sha256(receipt_path.read_bytes()).hexdigest()), 'review_script_sha256 = '+json.dumps(hashlib.sha256(Path(__file__).read_bytes()).hexdigest())]
(root/'PRODUCTION_PASS.toml').write_text('\n'.join(lines)+'\n')

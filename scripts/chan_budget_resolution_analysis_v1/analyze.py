"""Read-only, small-array analysis; execute on Unity after all three branches finish."""
import csv, json, sys, tomllib
from pathlib import Path
import numpy as np

root=Path('/work/pi_nicholas_pizzo_uri_edu/jupiter-ffr-breeze')
job=sys.argv[1]
result={'gpu_array_job':job,'cases':{},'interpretation_limit':'Short local budget sensitivity; no long-term LSV convergence claim. Spatial transfer can seed adjustment; compare full and late intervals.'}
for i,name in enumerate(('baseline','refined','smaller_dt')):
    d=root/f'data/chan-budget-resolution-v2-run-{job}_{i}/branch'
    meta=tomllib.loads((d/'metadata.toml').read_text())
    assert meta['status']=='COMPLETE', (name,meta['status'])
    with (d/'diagnostics.csv').open() as f: rows=list(csv.DictReader(f))
    cols={k:np.array([float(r[k]) for r in rows]) for k in rows[0]}
    t=cols['time']; elapsed=t-t[0]; residual=cols['residual_step']
    assert elapsed[-1]>=9.999
    rates={}
    for start in (0,1,2):
        mask=elapsed>=start-1e-7
        slope=np.polyfit(elapsed[mask],residual[mask],1)[0]
        rates[str(start)]={'residual_slope':float(slope),'slope_over_bottom_power':float(slope/.5)}
    sparse=np.abs(elapsed-np.rint(elapsed))<1e-7
    result['cases'][name]={
        'time_start':float(t[0]),'time_end':float(t[-1]),'grid':[meta['nx'],meta['ny'],meta['nz']],
        'initial_energy':float(cols['energy'][0]),'final_incremental_residual':float(residual[-1]),
        'residual_over_bottom_heat':float(residual[-1]/cols['bottom_heat'][-1]),
        'residual_rates_by_excluded_initial_time':rates,
        'heat_sampling_difference_final':float(cols['heat_sparse'][-1]-cols['heat_step'][-1]),
        'heat_sampling_difference_max_at_common_samples':float(np.max(np.abs(cols['heat_sparse'][sparse]-cols['heat_step'][sparse]))),
        'max_mach':float(cols['max_mach'].max()),'max_mass_error':float(cols['mass_error'].max()),
        'reported_dt_min':float(cols['dt'].min()),'reported_dt_max':float(cols['dt'].max()),
        'wall_seconds':meta['wall_seconds_actual']}
    if name=='baseline':
        parent=root/'data/chan-clean-pair-v1-run-2214011_0/parent/diagnostics.csv'
        with parent.open() as f: p=list(csv.DictReader(f))
        match=[r for r in p if abs(float(r['time'])-t[0])<1e-7]
        if match: result['baseline_initial_energy_minus_parent']=float(cols['energy'][0]-float(match[0]['energy_proxy']))
out=root/'runs/chan-budget-resolution-v2'/f'analysis-{job}.json'
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

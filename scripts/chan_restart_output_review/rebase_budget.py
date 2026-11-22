"""Unity-only derived correction; never overwrite raw output or simulation sources.
The manual output immediately after checkpoint restore preceded derived-state
initialization. Exclude that point/frame and rebase budgets on the first post-step
sample. This corrects reporting only; it does not repair or alter evolution.
"""
from pathlib import Path
import csv,json,sys
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((root/'diagnostics.csv').open())]
assert rows[0]['time']==263.06682088232816
rows=[r for r in rows if r['time']>263.06682088232816];assert rows and rows[0]['time']==263.5
first=rows[0];heat=0.;previous=None
for r in rows:
 if previous is not None:heat+=(r['time']-previous['time'])*.5*(r['external_power']+previous['external_power'])
 r['cumulative_external_heat_since_restart']=heat
 r['cumulative_energy_defect_since_restart']=r['total_energy']-first['total_energy']-heat
 r['energy_defect_over_restart_ke']=r['cumulative_energy_defect_since_restart']/first['kinetic_energy']
 if previous is None:
  r['energy_residual_over_input']=float('nan');r['centered_residual_over_input']=float('nan')
 previous=r
with (out/'diagnostics_from_first_poststep.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
record={'excluded_time':263.06682088232816,'excluded_movie_frame':'frame_0026307.jld2','valid_budget_start_time':263.5,'last_time':rows[-1]['time'],'reason':'Manual restart output preceded the derived-state update. Restored prognostic state is advanced normally by run initialization; CPU control comparison matched previous implementation exactly. Raw first output and cumulative offsets are invalid diagnostics, not a physical energy jump.','missing_budget_interval':[263.06682088232816,263.5],'corrected_final_energy_defect':rows[-1]['cumulative_energy_defect_since_restart'],'corrected_net_heat':heat,'raw_output_modified':False,'recompute_after_completion':True}
(out/'CORRECTION.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record))

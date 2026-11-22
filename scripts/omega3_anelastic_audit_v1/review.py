import sys,json,csv,tomllib
from pathlib import Path
import numpy as np,h5py
root=Path('/work/pi_nicholas_pizzo_uri_edu/jupiter-ffr-breeze');out=root/sys.argv[1]
a=root/sys.argv[2];c=root/'data/chan-eddies-v4-case4-2090026/parent'
assert tomllib.loads((a/'metadata.toml').read_text())['status']=='COMPLETE'
ma=tomllib.loads((a/'INITIAL_MAPPING.toml').read_text())
report={'initial_mapping':ma,'pilot':str(a.relative_to(root)),'compressible_control':str(c.relative_to(root)),'long_run_released':False}
rows=list(csv.DictReader((a/'diagnostics.csv').open()));r0,r1=rows[0],rows[-1]
report['pilot_summary']={k:float(r1[k]) for k in ['time','wall_seconds','max_mach','constraint_relative','theta_departure_max','thermodynamic_density_departure_max','kinetic_energy','urms','cooling_power']}
report['pilot_summary']['ke_ratio_final_initial']=float(r1['kinetic_energy'])/float(r0['kinetic_energy'])
report['pilot_summary']['maximum_mach']=max(float(r['max_mach']) for r in rows)
report['pilot_summary']['maximum_constraint_relative']=max(float(r['constraint_relative']) for r in rows)
late=[r for r in rows if float(r['time'])>=2.5]
rate=(float(late[-1]['wall_seconds'])-float(late[0]['wall_seconds']))/(float(late[-1]['time'])-float(late[0]['time']))
report['pilot_summary']['seconds_per_time_unit_late']=rate
report['pilot_summary']['estimated_hours_t5_to1000_at_pilot_rate']=rate*995/3600
cr=list(csv.DictReader((c/'diagnostics.csv').open()))
near=min(cr,key=lambda x:abs(float(x['time'])-5))
report['compressible_at_nearest_time']={k:float(near[k]) for k in ['time','max_mach','kinetic_energy','urms_conv','cooling_power']}
report['compressible_at_nearest_time']['ke_ratio_final_initial']=float(near['kinetic_energy'])/float(ma['compressible_initialization']['kinetic_energy'])
report['interpretation']='Validation only. Matched velocity and entropy do not imply matched finite-amplitude thermodynamic pressure/density; fixed-reference and compressible equations differ. No vortex-selection conclusion from t<=5.'
out.mkdir(parents=True,exist_ok=False);(out/'review.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

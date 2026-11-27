"""Scientific release criteria; run on Unity after B0, never on local fields."""
from pathlib import Path
import csv,json,sys,tomllib
import numpy as np,h5py
root=Path.cwd();src=Path(sys.argv[1]);out=root/'runs/cai2022-benchmark-v1'
meta=tomllib.loads((src/'metadata.toml').read_text());rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((src/'diagnostics.csv').open())]
late=[r for r in rows if r['time']>=2000];result={'status':'NOT_READY','source':str(src),'time_end':rows[-1]['time'],'case':'B0','criteria_definition':'Proposed benchmark acceptance tolerances, not published error bars.'}
if meta['status']!='COMPLETE' or rows[-1]['time']<3000 or len(late)<50:
 result['reason']='B0 thermal-calibration integration incomplete; no B7 production release.'
else:
 t=np.array([r['time'] for r in late]);u=np.array([r['urms_conv'] for r in late]);q=np.array([r['top_flux'] for r in late]);pr=np.array([r['Pr_mean'] for r in late]);re=np.array([r['Re_mean'] for r in late])
 trend=float(np.polyfit(t,u,1)[0]*1000/u.mean());flux_ratio=float(q.mean()/meta['bottom_flux']);flux=[]
 for file in sorted(src.glob('frame_*.jld2')):
  if int(file.stem.split('_')[1])/1000<2000:continue
  with h5py.File(file) as h:
   z=h['z'][:];flux.append(float(h['enthalpy_flux_profile'][:][np.argmin(abs(z-.5))]))
 conv_ratio=float(np.mean(flux)/meta['bottom_flux'])
 checks={'rms_within_25pct_of_0.066':.0495<=u.mean()<=.0825,'top_bottom_flux_balance_10pct':.9<=flux_ratio<=1.1,'rms_trend_less_10pct_per_thermal_time':abs(trend)<.1,'mid_layer_enthalpy_fraction_between_0.1_and_0.3':.1<=conv_ratio<=.3,'Pr_within_factor_1.5_of_0.396':.396/1.5<=pr.mean()<=.396*1.5,'Re_within_factor_1.5_of_4703':4703/1.5<=re.mean()<=4703*1.5,'energy_residual_below_0.1pct':max(abs(r['defect_over_bottom_heat']) for r in late)<.001}
 result.update(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,mean_rms=float(u.mean()),rms_trend_per_thermal_time=trend,mean_top_flux_over_bottom=flux_ratio,mean_mid_enthalpy_over_bottom=conv_ratio,mean_Pr=float(pr.mean()),mean_Re=float(re.mean()))
(out/'B0_calibration_assessment.json').write_text(json.dumps(result,indent=2)+'\n')
if result['status']=='PASS':(out/'B0_CALIBRATION_PASS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))

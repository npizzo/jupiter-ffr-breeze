"""Read-only Unity review of the short developed-flow resolution screen."""
import csv,json,sys
from pathlib import Path
import numpy as np
root=Path(sys.argv[1]);reports={}
for n in [128,256]:
 p=root/f'grid{n}'
 with open(p/'transfer.csv') as f: rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
 with open(p/'diagnostics.csv') as f: diag=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
 assert diag[-1]['time']>=12020-1e-8
 rr=[r for r in rows if r['time']>=12005]
 assert len(rr)>10
 rate=lambda key:np.array([r[key]/r['dt'] for r in rr])
 conv=rate('convective_stress_work');dis=rate('native_minus_centered_flux_work');res=rate('closure_residual')
 scale=max(float(np.mean(abs(conv))),1e-25)
 reports[str(n)]={'sampled_steps':len(rr),'mean_convective_work_rate':float(conv.mean()),'mean_abs_convective_work_rate':scale,'mean_native_discrepancy_rate':float(dis.mean()),'mean_abs_discrepancy_over_convective':float(np.mean(abs(dis))/scale),'max_closure_residual_over_convective':float(np.max(abs(res))/scale),'max_rhs_partition_relative':max(r['rhs_partition_relative'] for r in rr),'elapsed_wall':diag[-1]['wall_seconds'],'simulated_time':diag[-1]['time']-diag[0]['time'],'lambda_cells_final':diag[-1]['lambda_cells'],'final_urms':diag[-1]['urms_conv'],'cap_min_gradient_final':diag[-1]['cap_min_dtheta_dz']}
a,b=reports['128'],reports['256'];scale=max(a['mean_abs_convective_work_rate'],b['mean_abs_convective_work_rate'],1e-25)
difference=abs(a['mean_convective_work_rate']-b['mean_convective_work_rate'])/scale
checks={'native_closure':all(v['max_closure_residual_over_convective']<1e-4 for v in reports.values()),'resolved_transfer_reconstruction':all(v['mean_abs_discrepancy_over_convective']<.5 for v in reports.values()),'coarse_fine_transfer_screen':difference<.5}
result={'status':'PASS' if all(checks.values()) else 'REVIEW_REQUIRED_DO_NOT_RELEASE','checks':checks,'grids':reports,'coarse_fine_mean_transfer_difference_over_convective':difference,'interpretation':'Short numerical/transfer screen; not long-time convergence or a sign-selection result. Native WENO minus centered flux is measured separately.','cost_note':'Conservative extrapolation includes compilation and dense every-step observer overhead; inspect measured production cadence before release if this alone fails.'}
(root/'SCREEN_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
(root/'RESPONSE_RELEASE.toml').write_text('status = "'+result['status']+'"\n')
print(json.dumps(result,indent=2))

from pathlib import Path
import csv,json,tomllib,h5py,numpy as np
p=Path('data/beta-shear-physics-v1-run-2162603');result={'time_unit':'nondimensional','late_window':[75,100],'cases':{}}
for name in ['nojet','jet']:
 dirs=[p/(name+'_pilot'),p/(name+'_to100')];metrics={};diag={};extent={}
 for d in dirs:
  for row in csv.DictReader((d/'vortex_metrics.csv').open()):
   r={k:(v if k=='component' else float(v)) for k,v in row.items()};metrics[(r['time'],r['z'],r['component'])]=r
  for row in csv.DictReader((d/'diagnostics.csv').open()):diag[float(row['time'])]={k:float(v) for k,v in row.items()}
  for row in csv.DictReader((d/'hotspot_extent.csv').open()):extent[float(row['time'])]={k:float(v) for k,v in row.items()}
 case={'status':tomllib.loads((dirs[-1]/'metadata.toml').read_text())['status'],'last_diagnostics':diag[max(diag)],'maxima':{k:max(r[k] for r in extent.values()) for k in ['mach_max','outside_mach','bulk_mach','fraction_above04']},'late_metrics':{},'late_interior_metrics':{},'initial_metrics':[r for r in metrics.values() if r['time']==0]}
 for z in sorted({r['z'] for r in metrics.values()}):
  case['late_metrics'][str(z)]={}
  for comp in ['total','nonzonal']:
   rr=[r for r in metrics.values() if r['z']==z and r['component']==comp and 75<=r['time']<=100]
   case['late_metrics'][str(z)][comp]={k:float(np.mean([r[k] for r in rr])) for k in ['skewness','anticyclonic_enstrophy_share','rms_zeta_over_f']}
 interior={}
 for fp in sorted(dirs[-1].glob('frame_*.jld2')):
  with h5py.File(fp) as h:
   t=float(h['t'][()])
   if t<75:continue
   q=h['vorticity_slices'][:]/float(h['f'][()]);z=h['video_heights'][:];ly=float(h['ly'][()]);ny=q.shape[1];y=(np.arange(ny)+.5)*ly/ny-ly/2
   for k,height in enumerate(z):
    for comp in ['total','nonzonal']:
     a=q[k]
     if comp=='nonzonal':a=a-a.mean(axis=1,keepdims=True)
     a=a[np.abs(y)<3.5];mu=a.mean();var=np.mean((a-mu)**2)
     values={'skewness':float(np.mean((a-mu)**3)/max(var,1e-30)**1.5),'anticyclonic_enstrophy_share':float(np.sum(a*a*(a<0))/max(np.sum(a*a),1e-30))}
     interior.setdefault((str(float(height)),comp),[]).append(values)
 for (z,comp),rr in interior.items():case['late_interior_metrics'].setdefault(z,{})[comp]={k:float(np.mean([r[k] for r in rr])) for k in rr[0]}
 with h5py.File(dirs[-1]/'profiles_0010000.jld2') as h:
  z=h['z'][:];n2=h['mean_N2'][:];zf=.5*(z[1:]+z[:-1]);case['final_cap_mean_N2_min']=float(n2[zf>.82].min())
 # Record rather than assume continuous diagnostic budgets across restart.
 a=list(csv.DictReader((dirs[0]/'diagnostics.csv').open()))[-1];b=list(csv.DictReader((dirs[1]/'diagnostics.csv').open()))[0]
 case['restart_budget_difference']={k:float(b[k])-float(a[k]) for k in ['time','energy','cumulative_external_work','energy_defect']}
 result['cases'][name]=case
(p/'PHYSICS_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

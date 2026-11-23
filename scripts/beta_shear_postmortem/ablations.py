import csv,json,pathlib,tomllib
root=pathlib.Path('data/beta-shear-diag-v2-2136871')
out=pathlib.Path('data/beta-shear-v1-postmortem')
result={'source_job':2136871,'cases':{},'scope':'Short x-uniform cross-sections; diagnostic ablations, not production nor a converged nonlinear failure replay.'}
keep=['time','dt','min_temperature','max_mach','nonzonal_velocity_rms','jet_profile_rms_error','kinetic_energy','energy_defect_over_restart_ke']
for name in ['baseline','half_dt','no_entropy_diffusion','no_horizontal_rotation']:
 p=root/name
 if not (p/'diagnostics.csv').exists():
  result['cases'][name]={'status':'NO_OUTPUT'};continue
 rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader((p/'diagnostics.csv').open())]
 meta=tomllib.loads((p/'metadata.toml').read_text())
 result['cases'][name]={'status':meta['status'],'last':{k:rows[-1][k] for k in keep},'samples':[{k:r[k] for k in keep} for r in rows if r['time'] in [0,5,10,15,20,25,30]],'grid':[meta['nx'],meta['ny'],meta['nz']],'dx':meta['lx']/meta['nx'],'dy':meta['ly']/meta['ny'],'initial_rms':meta['initial_rms'],'perturbation':meta['perturbation'],'entropy_diffusion_factor':meta.get('entropy_diffusion_factor',1),'horizontal_rotation_factor':meta.get('horizontal_rotation_factor',1)}
(out/'ablation-summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)

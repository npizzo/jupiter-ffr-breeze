"""Read saved-output schemas on Unity; record what can be rendered without invented depths."""
from pathlib import Path
import json,tomllib,csv,h5py
rows=[]
for p in sorted(Path('data').glob('**/metadata.toml')):
 try:m=tomllib.loads(p.read_text())
 except Exception:continue
 frames=sorted(p.parent.glob('frame_*.jld2'))
 r={'source':str(p.parent),'metadata_status':m.get('status','unspecified'),'frame_count':len(frames),'three_saved_vorticity_planes':False,'saved_fields':[],'last_saved_day':None,'classification':'no_saved_frames'}
 if frames:
  try:
   with h5py.File(frames[0]) as f:
    r['saved_fields']=list(f.keys()); r['three_saved_vorticity_planes']='vorticity_slices' in f and f['vorticity_slices'].shape[0]==3
   with h5py.File(frames[-1]) as f:r['last_saved_day']=float(f['t'][()])/86400
   technical=any(w in str(p.parent) for w in ('restart','check_','smoke','cpu_smoke','quick_','continuous','first_segment','resumed'))
   r['classification']='technical_validation' if technical else ('three_depth_output' if r['three_saved_vorticity_planes'] else 'legacy_output_without_three_saved_vorticity_planes')
  except Exception as e:r['classification']='unreadable_frame';r['error']=str(e)
 rows.append(r)
out=Path('rendered/gallery-coverage-20261111');out.mkdir(exist_ok=True)
(out/'inventory.json').write_text(json.dumps({'compute_location':'Unity','note':'Inventory of saved files, not a convergence or physical-validation claim. Metadata RUNNING is not current scheduler status. Technical restart/smoke tests are catalogued separately from science movies.','cases':rows},indent=2)+'\n')
with (out/'inventory.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['source','metadata_status','frame_count','three_saved_vorticity_planes','last_saved_day','classification']);w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in rows)
print(json.dumps({'cases':len(rows),'legacy':[r['source'] for r in rows if r['classification']=='legacy_output_without_three_saved_vorticity_planes'],'unreadable':[r for r in rows if r['classification']=='unreadable_frame']}),flush=True)

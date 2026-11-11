"""Read only saved tracer planes on Unity; report late-time contrast and precision."""
from pathlib import Path
import h5py,numpy as np,json
root=Path.cwd();out=root/'rendered/tracer100-contrast';out.mkdir(exist_ok=True)
rows=[]
for mode in ['heated','control']:
 case=root/'data/cascade100-v1'/('cascade100_f_'+mode)
 for day in [0,10,25,50,75,90,100]:
  index=int(day*8);paths=list(case.glob(f'segment_*/frame_{index:05d}.jld2'));assert len(paths)==1,(day,paths)
  with h5py.File(paths[0]) as d:
   saved=d['tracer_slices'][:];values=saved.astype(np.float64)
   for level,c in enumerate(values):
    mean=float(c.mean());delta=c-mean
    rows.append({'case':mode,'day':day,'level':level,'mean':mean,'min':float(c.min()),'max':float(c.max()),'anomaly_abs_quantiles':np.quantile(np.abs(delta),[.5,.9,.99,1]).tolist(),'saved_dtype':str(saved.dtype),'unique_values':int(np.unique(saved[level]).size),'float32_spacing_near_mean':float(abs(np.spacing(np.float32(mean))))})
report={'compute_location':'Unity','status':'COMPLETE','rows':rows}
(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
for row in rows:
 if row['day']>=75:print(json.dumps(row),flush=True)

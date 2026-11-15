"""Reduce one temporary 3D export on Unity, then retire that scratch export."""
from pathlib import Path
import sys,json,hashlib,os
import numpy as np
from diagnostics import field_analysis
p=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=True,parents=True)
assert p.parent.name=='fields3d' and p.name.startswith('state_')
row,maps=field_analysis(p)
def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.generic):return clean(v.item())
 if isinstance(v,float) and not np.isfinite(v):return None
 return v
np.savez_compressed(out/(p.stem+'.npz'),**maps)
report=out/(p.stem+'.json');report.write_text(json.dumps(clean(row),allow_nan=False))
with np.load(out/(p.stem+'.npz')) as d:assert np.isfinite(d['Z']).all() and d['Z'].shape==maps['Z'].shape
p.unlink() # Deliberately temporary export; retain it on any reduction failure.
print('REDUCED_TEMPORARY_FIELD',p.name,flush=True)

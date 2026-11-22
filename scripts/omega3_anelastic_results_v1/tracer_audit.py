import json
from pathlib import Path
import h5py,numpy as np
out=Path('data/omega3-anelastic-results-v1-2120352/tracer_audit.json');assert not out.exists()
root=Path('data/omega3-anelastic-long-v1-run-2113514/parent');result={}
for t in [100,500,1000]:
 with h5py.File(root/f'frame_{t*100:07d}.jld2') as f:
  c=f['tracer_slices'][()].astype(float);z=f['video_heights'][()]
 result[t]=[dict(height=float(z[k]),minimum=float(c[k].min()),maximum=float(c[k].max()),std=float(c[k].std())) for k in range(3)]
out.write_text(json.dumps(result,indent=2)+'\n');print(out.read_text())

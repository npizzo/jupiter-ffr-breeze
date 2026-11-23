"""Read-only scientific release gate for this bounded heated pair, not production."""
import json,tomllib,sys
from pathlib import Path
import numpy as np,h5py
here=Path(__file__).resolve().parent;root=here.parents[1]
c=tomllib.loads((here/'candidate.toml').read_text());thermal=root/c['thermal_parent'];geometry=root/c['geometry_parent']
screens=tomllib.loads((geometry/'SCREENS.toml').read_text())
for sign in [-1,1]:
    s=screens[f"H{c['lz']}_U{sign*c['jet_speed']}"]
    assert s['accepted_static'] and s['Tmin']>.7
outcomes=tomllib.loads((thermal/'THERMAL_OUTCOMES.toml').read_text())
assert len(outcomes)==6
assert all(m['status']=='COMPLETE' and m['final_time']>=40-1e-8 for m in outcomes.values()),'Thermal controls did not all finish'
def profiles(path):
    d={}
    for p in path.glob('profiles_*.jld2'):
        with h5py.File(p) as h:d[round(float(h['t'][()]),8)]=h['mean_temperature'][:]
    return d
def effect(label):
    no=profiles(thermal/(label+'_nojet'));jet=profiles(thermal/(label+'_jet'))
    times=sorted(no.keys()&jet.keys());assert times[0]==0 and times[-1]>=40-1e-8
    initial=jet[0]-no[0]
    return max(float(np.max(abs(jet[t]-no[t]-initial)/no[0])) for t in times)
original=effect('original_entropy1');lilly=effect('lilly_entropy1')
assert lilly<original, 'Lilly did not reduce the matched jet-induced drift'
assert lilly<.05, 'Matched jet-induced thermal drift exceeds 5%'
r={'status':'PASS','original_matched_drift':original,'lilly_matched_drift':lilly,'production_release':False,
   'scope':'Release bounded t60 heated validation only; does not establish thermal equilibrium.'}
Path(sys.argv[1]).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)

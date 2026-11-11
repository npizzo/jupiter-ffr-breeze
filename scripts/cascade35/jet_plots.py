"""Matched dimensional zonal-flow histories; run only on Unity."""
from pathlib import Path
import sys,json,tomllib,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
cases=[Path(p) for p in sys.argv[1:3]];out=Path(sys.argv[3]);records=[]
for case in cases:
 files=[]
 for segment in sorted(case.glob('segment_*'),key=lambda p:int(p.name.split('_')[-1])):
  m=tomllib.loads((segment/'metadata.toml').read_text());assert m['status']=='COMPLETE'
  files.extend(sorted(segment.glob('frame_*.jld2')))
 times=[];jets=[];energies=[]
 for file in files:
  with h5py.File(file) as d:
   u=d['ubar'][:];v=d['vbar'][:];jet=u.mean(axis=1)
   times.append(float(d['t'][()]));jets.append(jet)
   energies.append(float(np.mean(jet**2)/max(np.mean(u**2+v**2),1e-30)))
 times=np.asarray(times);jets=np.asarray(jets);assert np.all(np.isfinite(jets))
 assert np.array_equal(times,np.arange(801)*10800)
 records.append({'latitude':m['latitude'],'times':times,'jets':jets,'fraction':np.asarray(energies),'ly':m['ly'],'ny':m['ny'],'source':str(case)})
assert np.array_equal(records[0]['times'],records[1]['times'])
lim=max(float(np.max(np.abs(r['jets']))) for r in records);lim=max(lim,1e-6)
fig,axes=plt.subplots(3,1,figsize=(12,10),layout='constrained')
for ax,r in zip(axes[:2],records):
 y=(np.arange(r['ny'])+.5)*r['ly']/r['ny']/1000-r['ly']/2000
 im=ax.pcolormesh(r['times']/86400,y,r['jets'].T,shading='nearest',cmap='RdBu_r',norm=AsinhNorm(linear_width=.5,vmin=-lim,vmax=lim))
 ax.set(xlabel='Earth day',ylabel='y [km]',title=f"{r['latitude']:g}°N beta-plane: density-weighted depth-mean zonal flow")
 fig.colorbar(im,ax=ax,label='Zonal velocity [m/s], shared fixed asinh scale')
 axes[2].plot(r['times']/86400,r['fraction'],label=f"{r['latitude']:g}°N")
axes[2].set(xlabel='Earth day',ylabel='Zonal-mean / total depth-mean KE',ylim=(0,1));axes[2].legend();axes[2].grid(alpha=.2)
fig.savefig(out/'jet_evolution.png',dpi=150);plt.close(fig)
np.savez_compressed(out/'jet_evolution.npz',times_seconds=records[0]['times'],latitude55_zonal_u=records[0]['jets'],latitude35_zonal_u=records[1]['jets'],latitude55_fraction=records[0]['fraction'],latitude35_fraction=records[1]['fraction'])
(out/'jet_evolution.json').write_text(json.dumps({'compute_location':'Unity','latitudes':[r['latitude'] for r in records],'source_cases':[r['source'] for r in records],'frames':801,'start_day':0,'end_day':100,'velocity_units':'m/s','color_limits':[-lim,lim],'asinh_width_m_s':.5,'energy_fraction_definition':'mean_y(mean_x(ubar)^2)/mean_xy(ubar^2+vbar^2); includes any uniform zonal flow; density-weighted full-depth means; not the projected spectral metric','note':'Persistent bands and growth of zonal energy are diagnostics of organization, not proof of equilibrated jets or Ferrel cells.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print('Verified jet histories at both latitudes',flush=True)

"""Re-render the whole-run mean with an unclipped title; run only on Unity."""
from pathlib import Path
import h5py,numpy as np,json,tomllib,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
source=Path('scripts/ferrel35_day300_audit.py').read_text()
# Reuse exactly the audit's visualization filter, without rerunning its analysis.
ns={'np':np};exec(source[source.index('def gaussian_filter'):source.index('roots=')],ns)
filter_=ns['gaussian_filter'];out=Path('rendered/ferrel35-day300-audit-v1')
summary=json.loads((out/'summary.json').read_text());means={};count=0
for seg in range(1,31):
 root=Path('data/cascade35-v1/cascade35_heated' if seg<=10 else 'data/cascade35-300-v1/cascade35_heated')/f'segment_{seg}'
 m=tomllib.loads((root/'metadata.toml').read_text())
 for p in sorted(root.glob('frame_*.jld2')):
  with h5py.File(p) as d:
   t=float(d['t'][()]);assert t==count*10800
   weight=.5 if count in (0,2400) else 1.
   for key in ['mass_streamfunction','zonal_u']:
    a=d[key][:]
    if key not in means:means[key]=np.zeros_like(a,dtype=float)
    means[key]+=weight*a/2400
  count+=1
assert count==2401
assert np.isclose(abs(means['mass_streamfunction']).max(),summary['windows'][-1]['psi_max_abs_kg_m_s'],rtol=1e-12)
ny,nz=m['ny'],m['nz'];dy=m['ly']/ny;dz=m['height']/nz
y=(np.arange(ny)+.5)*dy-m['ly']/2;z=(np.arange(nz)+.5)*dz
psi=filter_(means['mass_streamfunction'],(5000/dz,200000/dy))/1e5
u=filter_(means['zonal_u'],(5000/dz,200000/dy))
fig,axs=plt.subplots(2,1,figsize=(7,8),layout='constrained')
lim=abs(psi).max();im=axs[0].pcolormesh((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,psi,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
axs[0].contour((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,psi,levels=np.linspace(-lim,lim,13),colors='k',linewidths=.35)
axs[0].set(title='Days 0–300: mass streamfunction');fig.colorbar(im,ax=axs[0],label='Psi / 1e5 [kg m^-1 s^-1]')
ul=abs(u).max();im=axs[1].pcolormesh(y/1000,z/1000,u,cmap='RdBu_r',vmin=-ul,vmax=ul,shading='auto');axs[1].set(title='Days 0–300: zonal velocity');fig.colorbar(im,ax=axs[1],label='u [m/s]')
for ax in axs:ax.set(xlabel='y north [km]',ylabel='Height [km]');ax.axhline(80,color='k',ls='--',lw=.6)
fig.suptitle('35°N heated beta-plane | whole-run mean\nIncludes initial adjustment; not a stationary mean\nSmoothing: 200 km in y, 5 km in z; dashed: initial cap base',fontsize=12)
fig.savefig(out/'overturning_whole_run.png',dpi=130);plt.close(fig)
(out/'whole_run_plot_provenance.json').write_text(json.dumps({'status':'PASS','compute_location':'Unity','frames':count,'same_trapezoidal_mean':True,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
print('Whole-run figure repaired; same 2401-frame trapezoidal mean verified.')

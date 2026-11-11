"""Save 0-300 day zonal-flow histories and block statistics on Unity."""
from pathlib import Path
import sys,json,tomllib,csv
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
case,out=map(Path,sys.argv[1:3]);files=[]
for p in sorted(case.glob('segment_*'),key=lambda p:int(p.name.split('_')[-1])):
 m=tomllib.loads((p/'metadata.toml').read_text());assert m['status']=='COMPLETE'
 files+=sorted(p.glob('frame_*.jld2'))
times=[];jets=[];fractions=[];psi_max=[];mass_residual=[];native_flux=[]
for p in files:
 with h5py.File(p) as d:
  u=d['ubar'][:];v=d['vbar'][:];jet=u.mean(axis=1)
  times.append(float(d['t'][()]));jets.append(jet)
  fractions.append(float(np.mean(jet**2)/max(np.mean(u**2+v**2),1e-30)))
  psi_max.append(float(np.max(np.abs(d['mass_streamfunction'][:]))))
  mass_residual.append(float(d['mass_flux_reconstruction_max'][()]))
  native_flux.append(max(float(np.max(np.abs(d['native_rhow'][:]))),1e-30))
times=np.asarray(times);jets=np.asarray(jets);fractions=np.asarray(fractions)
assert np.array_equal(times,np.arange(2401)*10800) and np.isfinite(jets).all()
days=times/86400;y=(np.arange(m['ny'])+.5)*m['ly']/m['ny']/1000-m['ly']/2000
limit=max(float(np.max(np.abs(jets))),1e-6)
fig,axes=plt.subplots(3,1,figsize=(12,10),layout='constrained')
im=axes[0].pcolormesh(days,y,jets.T,shading='nearest',cmap='RdBu_r',norm=AsinhNorm(linear_width=.5,vmin=-limit,vmax=limit))
axes[0].set(xlabel='Earth day',ylabel='y [km]',title='35°N heated beta-plane: density-weighted depth-mean zonal flow')
fig.colorbar(im,ax=axes[0],label='Zonal velocity [m/s], fixed scale')
axes[1].plot(days,fractions);axes[1].set(xlabel='Earth day',ylabel='Zonal-mean / depth-mean KE',ylim=(0,1));axes[1].grid(alpha=.2)
blocks=[]
for first,last in [(100,150),(150,200),(200,250),(250,300)]:
 mask=(days>first)&(days<=last);profile=jets[mask].mean(axis=0)
 axes[2].plot(y,profile,label=f'Days {first}–{last}')
 blocks.append({'start_day_exclusive':first,'end_day_inclusive':last,'frames':int(mask.sum()),'mean_zonal_ke_fraction':float(fractions[mask].mean()),'mean_profile_rms_m_s':float(np.sqrt(np.mean(profile**2))),'mean_profile_m_s':profile.tolist()})
axes[2].set(xlabel='y [km]',ylabel='Mean zonal velocity [m/s]');axes[2].legend();axes[2].grid(alpha=.2)
fig.savefig(out/'jet_evolution_0_300days.png',dpi=150);plt.close(fig)
np.savez_compressed(out/'jet_evolution_0_300days.npz',times_seconds=times,y_km=y,zonal_u_m_s=jets,zonal_ke_fraction=fractions,mass_streamfunction_max=psi_max,mass_flux_reconstruction_max=mass_residual,native_vertical_flux_max=native_flux)
report={'status':'DIAGNOSTICS_COMPLETE_REVIEW_REQUIRED','compute_location':'Unity','frames':2401,'blocks':blocks,'late_profile_change_rms_m_s':float(np.sqrt(np.mean((np.asarray(blocks[-1]['mean_profile_m_s'])-np.asarray(blocks[-2]['mean_profile_m_s']))**2))),'zonal_energy_fraction_definition':'mean_y(mean_x(ubar)^2)/mean_xy(ubar^2+vbar^2), includes uniform zonal flow; not projected spectral metric','maximum_mass_reconstruction_error':max(mass_residual),'continuation_to_500':'Not submitted. Review physical/numerical checks and late-time jet drift before deciding.','limitations':'Block differences are descriptive, not independent-sample confidence intervals. Jets or overturning images do not demonstrate equilibrium, a resolved inverse cascade, or Ferrel cells.'}
(out/'day300-review.json').write_text(json.dumps(report,indent=2)+'\n')
print('Jet histories and review statistics complete',flush=True)

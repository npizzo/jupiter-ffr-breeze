"""Run on Unity only; audit the native-face decomposition and thermal budget."""
from pathlib import Path
import json,sys,tomllib,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
meta=tomllib.loads((root/'metadata.toml').read_text());assert meta['status']=='COMPLETE'
p=sorted(root.glob('frame_*.jld2'))[-1]
with h5py.File(p) as d:
 t0=float(d['face_origin_time'][()]);t1=float(d['face_end_time'][()]);dt=t1-t0
 assert t0==300*86400 and t1==301*86400 and int(d['face_steps'][()])==meta['native_budget_checks']
 rho=d['density'][:];y=d['y'][:];z=d['z'][:];pressure=d['pressure'][:]
 mom=d['face_integrated_momentum'][:]/dt/rho[None,:,None]
 native=d['native_integrated_parts'][:]/dt/rho[None,:,None]
 heat=d['face_integrated_heat'][:]/dt
 ht=d['face_integrated_heat_transport'][:]/dt
 state_rate=(d['face_heat_current'][:]-d['face_heat_origin'][:])/dt
 residual=d['face_heat_residual'][:]/dt
 checks={k:float(d[k][()]) for k in ['face_max_momentum_flux_error','face_max_heat_partition_error','face_max_heat_increment_error']}
 mask=(abs(y)[None,:]<=2400000)&(z[:,None]>=10000)&(z[:,None]<=75000)
 weights=np.broadcast_to(rho[:,None],mask.shape)[mask].copy();weights/=weights.sum()
 rms=lambda a:float(np.sqrt(np.sum(weights*a[mask]**2)))
 cosine=lambda a,b:float(np.sum(weights*a[mask]*b[mask])/np.sqrt(np.sum(weights*a[mask]**2)*np.sum(weights*b[mask]**2)))
 names=['mean_y','mean_z','eddy_y','eddy_z','reconstruction_y','reconstruction_z']
 errors={'momentum_face_sum_vs_native_rms_m_s2':rms(mom.sum(0)-native[0]),'heat_face_sum_vs_native_max_rho_theta_s':float(abs(ht.sum(0)-heat[0]).max()),'heat_increment_residual_max_rho_theta_s':float(abs(residual).max())}
 assert errors['momentum_face_sum_vs_native_rms_m_s2']<1e-12
 assert errors['heat_increment_residual_max_rho_theta_s']<1e-12
 eddy=mom[2]+mom[3];mean=mom[0]+mom[1];remainder=mom[4]+mom[5]
 exner=(pressure/1e5)**((8.314462618/meta['molar_mass'])/meta['heat_capacity'])
 dz=meta['height']/meta['nz'];energy=meta['heat_capacity']*exner[:,None]
 column=lambda a:float(dz*np.mean(np.sum(energy*a,axis=0)))
 hnames=['advection','diffusion','stable_cap','thermal_feedback','bottom_flux']
 thermal={k:column(a) for k,a in zip(hnames,heat)}
 thermal.update(storage=column(state_rate),residual=column(residual))
 report={'status':'PASS','compute_location':'Unity','start_day':300,'end_day':301,'native_steps':int(d['face_steps'][()]),'mask':{'abs_y_max_km':2400,'z_km':[10,75]},'native_acceleration_rms_m_s2':dict(zip(['advection','coriolis','diffusion','external_forcing'],map(rms,native))),'face_acceleration_rms_m_s2':dict(zip(names,map(rms,mom))),'combined_acceleration_rms_m_s2':{'mean':rms(mean),'eddy':rms(eddy),'reconstruction':rms(remainder)},'reconstruction_to_eddy_rms':rms(remainder)/rms(eddy),'eddy_cosine_against_negative_coriolis':cosine(eddy,-native[1]),'native_checks':checks,'closure_checks':errors,'column_enthalpy_budget_W_m2':thermal,'limitations':['One mature-state day tests diagnostic accounting; it is not a statistical Ferrel-cell identification.','Mean and eddy fluxes use linear interpolation onto the actual transport faces; reconstruction is native WENO minus that specified baseline, not necessarily dissipation.','All terms are integrated at every SSPRK3 stage. The independent momentum and potential-temperature state increments check closure.','Column enthalpy uses cp times prescribed reference Exner times rho-theta; advection and diffusion may change its volume integral. It is not a full kinetic plus internal plus gravitational energy budget.','Wall, cap and resolution sensitivity remain separate tests.'],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'review.json').write_text(json.dumps(report,indent=2)+'\n')
fig,axs=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
for ax,a,title in zip(axs.flat,[mean,eddy,remainder,native[1],native[2],mom.sum(0)-native[0]],['Mean transport','Eddy transport','WENO minus linear product','Coriolis','Explicit diffusion','Face sum minus native advection']):
 lim=max(float(np.max(abs(a[mask]))),1e-18);im=ax.pcolormesh(y/1e6,z/1000,a,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto');fig.colorbar(im,ax=ax,label='m s⁻²');ax.set(title=title,xlabel='Northward distance (1000 km)',ylabel='Height (km)',ylim=(10,75),xlim=(-2.4,2.4))
fig.suptitle('35°N mature-state diagnostic | Earth days 300–301 | native-stage averages\nIndependent color ranges; numerical accounting test, not a cell identification')
fig.savefig(out/'momentum_budget.png',dpi=130);plt.close(fig)
print(json.dumps(report,indent=2),flush=True)

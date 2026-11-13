"""Read-only 200-day mean-circulation audit; all computations on Unity."""
from pathlib import Path
import json,tomllib,hashlib
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def gaussian_filter(a,sigma,mode='reflect'):
    # Separable Gaussian visualization filter with reflected cell values.
    a=np.asarray(a,dtype=float)
    for axis,width in enumerate(sigma):
        radius=int(np.ceil(4*width));x=np.arange(-radius,radius+1)
        kernel=np.exp(-.5*(x/width)**2);kernel/=kernel.sum()
        a=np.apply_along_axis(lambda v:np.convolve(np.pad(v,radius,mode='symmetric'),kernel,mode='valid'),axis,a)
    return a
roots=[Path('data/cascade35-v1/cascade35_heated'), Path('data/cascade35-300-v1/cascade35_heated')];out=Path('rendered/ferrel35-day200-audit-v2');out.mkdir(exist_ok=False)
metadata=[];physical_metadata=[]
files={}
for seg in range(1,21):
 p=roots[int(seg>10)]/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text());assert m['status']=='COMPLETE'
 physical_metadata.append(m)
 metadata.append({'segment':seg,'path':str(p),'sha256':hashlib.sha256((p/'metadata.toml').read_bytes()).hexdigest(),'configuration_sha256':m['configuration_sha256']})
 for f in sorted(p.glob('frame_*.jld2')):files[int(f.stem.split('_')[1])]=f
assert sorted(files)==list(range(1601))
# The continuation config has a different stop time and hash; check physical settings directly.
physical_keys=['nx','ny','nz','lx','ly','height','f','beta','latitude','gravity','theta0','heat_flux','heat_capacity','molar_mass','base_pressure','cap_n2','cap_height','cap_min_fraction','cap_core_bottom','cap_core_top','cooling_time','sponge_time','viscosity','thermal_diffusivity','tracer_diffusivity','cap_diffusivity','mixing_model','balance_heat','cooling_model','output_seconds','script_sha256','controller_sha256']
for mm in physical_metadata:
 for k in physical_keys:assert mm[k]==physical_metadata[0][k],k
keys=['zonal_u','zonal_v','zonal_w','zonal_theta','eddy_uv','eddy_uw','mass_streamfunction','native_rhov','native_rhow']
ny,nz=m['ny'],m['nz'];dy=m['ly']/ny;dz=m['height']/nz
y=(np.arange(ny)+.5)*dy-m['ly']/2;z=(np.arange(nz)+.5)*dz
fcor=m['f']+m['beta']*y
fields={k:[] for k in keys};times=[];qa={k:[] for k in ['mass_divergence_max','mass_flux_reconstruction_max','top_streamfunction_max']}
jet_fractions=[]
for i in range(400,1601):
 with h5py.File(files[i]) as d:
  times.append(float(d['t'][()]))
  for k in keys:fields[k].append(d[k][:]);assert np.isfinite(fields[k][-1]).all(),(i,k)
  for k in qa:qa[k].append(float(d[k][()]))
  rho=d['density'][:]
  ubar,vbar=d['ubar'][:],d['vbar'][:]
  jet_fractions.append(float(np.mean(ubar.mean(axis=1)**2)/np.mean(ubar**2+vbar**2)))
fields={k:np.array(v) for k,v in fields.items()};times=np.array(times)
assert np.array_equal(times,np.arange(400,1601)*10800)
mask=(z[:,None]>=10000)&(z[:,None]<=75000)&(abs(y[None,:])<=.4*m['ly'])
weights=np.broadcast_to(rho[:,None],(nz,ny))[mask];weights=weights/weights.sum()
def rms(a):return float(np.sqrt(np.sum(weights*a[mask]**2)))
def cosine(a,b):return float(np.sum(weights*a[mask]*b[mask])/max(rms(a)*rms(b),1e-30))
def mean(a):return np.trapezoid(a,dx=10800,axis=0)/(10800*(len(a)-1))
def average(a):return {k:mean(v[a]) for k,v in fields.items()}
windows=[]
for start in [50,100,150]:
 for a,b in [(start,start+25),(start+25,start+50),(start,start+50)]:
  windows.append((f'{a}–{b}',slice((a-50)*8,(b-50)*8+1)))
means=[average(s) for _,s in windows]
reports=[]
for (label,sl),avg in zip(windows,means):
 u,v,w=[fields[k][sl] for k in ['zonal_u','zonal_v','zonal_w']]
 dt=(len(u)-1)*10800; tendency=(u[-1]-u[0])/dt
 eh=-np.gradient(avg['eddy_uv'],dy,axis=1,edge_order=2)
 ev=-np.gradient(rho[:,None]*avg['eddy_uw'],dz,axis=0,edge_order=2)/rho[:,None]
 coriolis=fcor[None,:]*avg['zonal_v']
 adv=mean(-v*np.gradient(u,dy,axis=2,edge_order=2)-w*np.gradient(u,dz,axis=1,edge_order=2))
 residual=tendency-coriolis-eh-ev-adv
 coherent_v=rms(avg['zonal_v']);coherent_w=rms(avg['zonal_w'])
 reports.append({'window_earth_days':label,'mean_v_rms_m_s':coherent_v,'mean_w_rms_m_s':coherent_w,'v_persistence_ratio':coherent_v/max(float(np.sqrt(mean(np.array([rms(a)**2 for a in v])))),1e-30),'w_persistence_ratio':coherent_w/max(float(np.sqrt(mean(np.array([rms(a)**2 for a in w])))),1e-30),'mean_u_rms_m_s':rms(avg['zonal_u']),'psi_max_abs_kg_m_s':float(abs(avg['mass_streamfunction']).max()),'acceleration_rms_m_s2':{k:rms(a) for k,a in [('du_dt',tendency),('f_v',coriolis),('eddy_horizontal',eh),('eddy_vertical',ev),('mean_advection',adv),('unclosed_residual',residual)]},'eddy_vs_minus_fv_cosine':cosine(eh+ev,-coriolis),'horizontal_eddy_vs_minus_fv_cosine':cosine(eh,-coriolis),'vertical_eddy_vs_minus_fv_cosine':cosine(ev,-coriolis),'residual_over_fv_rms':rms(residual)/max(rms(coriolis),1e-30)})
# Smooth only for visual inspection and a separately labeled robustness diagnostic.
# Physical widths 200 km meridional and 5 km vertical, not tuned to the output.
def smooth(a):return gaussian_filter(a,sigma=(5000/dz,200000/dy),mode='reflect')
for report,avg in zip(reports,means):
 eh=-np.gradient(avg['eddy_uv'],dy,axis=1,edge_order=2);ev=-np.gradient(rho[:,None]*avg['eddy_uw'],dz,axis=0,edge_order=2)/rho[:,None]
 report['smoothed_eddy_vs_minus_fv_cosine']=cosine(smooth(eh+ev),smooth(-fcor[None,:]*avg['zonal_v']))
psi=[a['mass_streamfunction'] for a in means]
psimask=(np.arange(nz+1)[:,None]*dz>=10000)&(np.arange(nz+1)[:,None]*dz<=75000)&(abs((np.arange(ny+1)[None,:]*dy-m['ly']/2))<=.4*m['ly'])
def cosflat(a,b):return float(np.sum(a[psimask]*b[psimask])/max(np.linalg.norm(a[psimask])*np.linalg.norm(b[psimask]),1e-30))
summary={'status':'COMPLETE','compute_location':'Unity','sources':[str(r) for r in roots],'segment_metadata':metadata,'end_day':200,'frames_analyzed':1201,'windows':reports,'half_window_streamfunction_cosine':cosflat(psi[6],psi[7]),'half_window_smoothed_streamfunction_cosine':cosflat(smooth(psi[6]),smooth(psi[7])),'half_window_v_cosine':cosine(means[6]['zonal_v'],means[7]['zonal_v']),'half_window_w_cosine':cosine(means[6]['zonal_w'],means[7]['zonal_w']),'jet_ke_fraction_last50_mean':float(np.mean(jet_fractions[-401:])),'jet_ke_fraction_day200':jet_fractions[-1],'quality_maxima':{k:max(v) for k,v in qa.items()},'metric_region':'|y| <= 2400 km, 10 <= z <= 75 km; density weighting for centered field norms','streamfunction_definition':'rho*v=-dPsi/dz; rho*w=dPsi/dy; kg/(m s); Eulerian zonal mean','limitations':['Centered resolved covariance budget, not exact discrete momentum closure. Residual includes viscosity, numerical stresses and diagnostic errors.','No test of residual/Lagrangian circulation or thermodynamic indirectness.','Time-window persistence is a descriptive check, not a statistical significance test.','Stable cap, sponge, free-slip rigid walls and imposed thermal feedback can affect circulation.'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
summary['half_window_period_days']=[150,200]
summary['fifty_day_comparison']=[{'days':windows[i][0],'streamfunction_cosine_vs_previous':cosflat(psi[i-3],psi[i]) if i>2 else None,'zonal_u_cosine_vs_previous':cosine(means[i-3]['zonal_u'],means[i]['zonal_u']) if i>2 else None} for i in [2,5,8]]
summary['reproduction_check']={}
old=json.loads(Path('rendered/ferrel35-day100-audit/summary.json').read_text())
for i in range(3):
 for k in ['residual_over_fv_rms','eddy_vs_minus_fv_cosine','psi_max_abs_kg_m_s']:
  assert np.isclose(reports[i][k],old['windows'][i][k],rtol=1e-12,atol=1e-14),(i,k)
summary['reproduction_check']['earlier_50_100_day_audit']='PASS'
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
def make_figure(indices,filename):
    fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    ps=[smooth(psi[i]) for i in indices];lim=max(abs(p).max() for p in ps)/1e5
    ul=max(abs(smooth(a['zonal_u'])).max() for a in [means[i] for i in indices])
    for c,((label,_),p,a) in enumerate(zip([windows[i] for i in indices],ps,[means[i] for i in indices])):
     im=axs[0,c].pcolormesh((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
     axs[0,c].contour((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,levels=np.linspace(-lim,lim,13),colors='k',linewidths=.35)
     axs[0,c].set(title=f'Days {label}: mass streamfunction',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[0,c],label='Psi / 1e5 [kg m^-1 s^-1]')
     im=axs[1,c].pcolormesh(y/1000,z/1000,smooth(a['zonal_u']),cmap='RdBu_r',vmin=-ul,vmax=ul,shading='auto');axs[1,c].set(title=f'Days {label}: zonal velocity',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[1,c],label='u [m/s]')
     for ax in axs[:,c]:ax.axhline(80,color='k',ls='--',lw=.6)
    fig.suptitle('35N heated beta-plane | time-mean overturning audit\nVisualization smoothing: 200 km in y, 5 km in z; dashed line = initial cap base')
    fig.savefig(out/filename,dpi=130);plt.close(fig)
make_figure([6,7,8],'overturning_audit.png')
make_figure([2,5,8],'overturning_evolution.png')
print(json.dumps(summary,indent=2),flush=True)

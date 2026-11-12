"""Read-only 100-day mean-circulation audit; all computations on Unity."""
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
root=Path('data/cascade35-v1/cascade35_heated');out=Path('rendered/ferrel35-day100-audit');out.mkdir(exist_ok=False)
files={}
for seg in range(1,11):
 p=root/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text());assert m['status']=='COMPLETE'
 for f in sorted(p.glob('frame_*.jld2')):files[int(f.stem.split('_')[1])]=f
assert sorted(files)==list(range(801))
keys=['zonal_u','zonal_v','zonal_w','zonal_theta','eddy_uv','eddy_uw','mass_streamfunction','native_rhov','native_rhow']
ny,nz=m['ny'],m['nz'];dy=m['ly']/ny;dz=m['height']/nz
y=(np.arange(ny)+.5)*dy-m['ly']/2;z=(np.arange(nz)+.5)*dz
fcor=m['f']+m['beta']*y
fields={k:[] for k in keys};times=[];qa={k:[] for k in ['mass_divergence_max','mass_flux_reconstruction_max','top_streamfunction_max']}
jet_fractions=[]
for i in range(400,801):
 with h5py.File(files[i]) as d:
  times.append(float(d['t'][()]))
  for k in keys:fields[k].append(d[k][:])
  for k in qa:qa[k].append(float(d[k][()]))
  rho=d['density'][:]
  ubar,vbar=d['ubar'][:],d['vbar'][:]
  jet_fractions.append(float(np.mean(ubar.mean(axis=1)**2)/np.mean(ubar**2+vbar**2)))
fields={k:np.array(v) for k,v in fields.items()};times=np.array(times)
assert np.array_equal(times,np.arange(400,801)*10800)
mask=(z[:,None]>=10000)&(z[:,None]<=75000)&(abs(y[None,:])<=.4*m['ly'])
weights=np.broadcast_to(rho[:,None],(nz,ny))[mask];weights=weights/weights.sum()
def rms(a):return float(np.sqrt(np.sum(weights*a[mask]**2)))
def cosine(a,b):return float(np.sum(weights*a[mask]*b[mask])/max(rms(a)*rms(b),1e-30))
def mean(a):return np.trapezoid(a,dx=10800,axis=0)/(10800*(len(a)-1))
def average(a):return {k:mean(v[a]) for k,v in fields.items()}
windows=[('50–75',slice(0,201)),('75–100',slice(200,401)),('50–100',slice(0,401))]
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
summary={'status':'COMPLETE','compute_location':'Unity','source':str(root),'end_day':100,'frames_analyzed':401,'windows':reports,'half_window_streamfunction_cosine':cosflat(psi[0],psi[1]),'half_window_smoothed_streamfunction_cosine':cosflat(smooth(psi[0]),smooth(psi[1])),'half_window_v_cosine':cosine(means[0]['zonal_v'],means[1]['zonal_v']),'half_window_w_cosine':cosine(means[0]['zonal_w'],means[1]['zonal_w']),'jet_ke_fraction_last50_mean':float(np.mean(jet_fractions)),'jet_ke_fraction_day100':jet_fractions[-1],'quality_maxima':{k:max(v) for k,v in qa.items()},'metric_region':'|y| <= 2400 km, 10 <= z <= 75 km; density weighting for centered field norms','streamfunction_definition':'rho*v=-dPsi/dz; rho*w=dPsi/dy; kg/(m s); Eulerian zonal mean','limitations':['Centered resolved covariance budget, not exact discrete momentum closure. Residual includes viscosity, numerical stresses and diagnostic errors.','No test of residual/Lagrangian circulation or thermodynamic indirectness.','Time-window persistence is a descriptive check, not a statistical significance test.','Stable cap, sponge, free-slip rigid walls and imposed thermal feedback can affect circulation.'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
fig,axs=plt.subplots(2,3,figsize=(15,8),layout='constrained')
ps=[smooth(p) for p in psi];lim=max(abs(p).max() for p in ps)/1e5
ul=max(abs(smooth(a['zonal_u'])).max() for a in means)
for c,((label,_),p,a) in enumerate(zip(windows,ps,means)):
 im=axs[0,c].pcolormesh((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
 axs[0,c].contour((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,levels=np.linspace(-lim,lim,13),colors='k',linewidths=.35)
 axs[0,c].set(title=f'Days {label}: mass streamfunction',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[0,c],label='Psi / 1e5 [kg m^-1 s^-1]')
 im=axs[1,c].pcolormesh(y/1000,z/1000,smooth(a['zonal_u']),cmap='RdBu_r',vmin=-ul,vmax=ul,shading='auto');axs[1,c].set(title=f'Days {label}: zonal velocity',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[1,c],label='u [m/s]')
 for ax in axs[:,c]:ax.axhline(80,color='k',ls='--',lw=.6)
fig.suptitle('35N heated beta-plane | time-mean overturning audit\nVisualization smoothing: 200 km in y, 5 km in z; dashed line = initial cap base')
fig.savefig(out/'overturning_audit.png',dpi=130);plt.close(fig)
print(json.dumps(summary,indent=2),flush=True)

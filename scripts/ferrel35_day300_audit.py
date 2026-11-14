"""Read-only 0–300-day mean-circulation audit; all computations on Unity."""
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
roots=[Path('data/cascade35-v1/cascade35_heated'), Path('data/cascade35-300-v1/cascade35_heated')];out=Path('rendered/ferrel35-day300-audit-v1');out.mkdir(exist_ok=False)
metadata=[];physical_metadata=[]
files={}
for seg in range(1,31):
 p=roots[int(seg>10)]/f'segment_{seg}';m=tomllib.loads((p/'metadata.toml').read_text());assert m['status']=='COMPLETE'
 physical_metadata.append(m)
 metadata.append({'segment':seg,'path':str(p),'sha256':hashlib.sha256((p/'metadata.toml').read_bytes()).hexdigest(),'configuration_sha256':m['configuration_sha256']})
 assert m['simulated_seconds']==seg*10*86400
 for f in sorted(p.glob('frame_*.jld2')):
  number=int(f.stem.split('_')[1]);assert number not in files,('duplicate_saved_time',number)
  files[number]=f
assert sorted(files)==list(range(2401))
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
for i in range(0,2401):
 with h5py.File(files[i]) as d:
  times.append(float(d['t'][()]))
  for k in keys:fields[k].append(d[k][:]);assert np.isfinite(fields[k][-1]).all(),(i,k)
  for k in qa:qa[k].append(float(d[k][()]))
  frame_rho=d['density'][:];assert np.isfinite(frame_rho).all() and (frame_rho>0).all()
  if i==0:rho=frame_rho
  else:assert np.array_equal(rho,frame_rho),('reference_density_changed',i)
  ubar,vbar=d['ubar'][:],d['vbar'][:]
  jet_fractions.append(float(np.mean(ubar.mean(axis=1)**2)/np.mean(ubar**2+vbar**2)))
fields={k:np.array(v) for k,v in fields.items()};times=np.array(times)
assert np.array_equal(times,np.arange(0,2401)*10800)
mask=(z[:,None]>=10000)&(z[:,None]<=75000)&(abs(y[None,:])<=2400000)
weights=np.broadcast_to(rho[:,None],(nz,ny))[mask];weights=weights/weights.sum()
def rms(a):return float(np.sqrt(np.sum(weights*a[mask]**2)))
def cosine(a,b):return float(np.sum(weights*a[mask]*b[mask])/max(rms(a)*rms(b),1e-30))
def mean(a):return np.trapezoid(a,dx=10800,axis=0)/(10800*(len(a)-1))
def average(a):return {k:mean(v[a]) for k,v in fields.items()}
windows=[]
for start in [0,50,100,150,200,250]:
 for a,b in [(start,start+25),(start+25,start+50),(start,start+50)]:
  windows.append((f'{a}–{b}',slice(a*8,b*8+1)))
windows.append(('0–300',slice(0,2401)))
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
psimask=(np.arange(nz+1)[:,None]*dz>=10000)&(np.arange(nz+1)[:,None]*dz<=75000)&(abs((np.arange(ny+1)[None,:]*dy-m['ly']/2))<=2400000)
def cosflat(a,b):return float(np.sum(a[psimask]*b[psimask])/max(np.linalg.norm(a[psimask])*np.linalg.norm(b[psimask]),1e-30))
summary={'status':'COMPLETE','compute_location':'Unity','sources':[str(r) for r in roots],'segment_metadata':metadata,'start_day':0,'end_day':300,'frames_analyzed':2401,'windows':reports,'half_window_streamfunction_cosine':cosflat(psi[15],psi[16]),'half_window_smoothed_streamfunction_cosine':cosflat(smooth(psi[15]),smooth(psi[16])),'half_window_v_cosine':cosine(means[15]['zonal_v'],means[16]['zonal_v']),'half_window_w_cosine':cosine(means[15]['zonal_w'],means[16]['zonal_w']),'jet_ke_fraction_last50_mean':float(np.mean(jet_fractions[-401:])),'jet_ke_fraction_day300':jet_fractions[-1],'quality_maxima':{k:max(v) for k,v in qa.items()},'metric_region':'|y| <= 2400 km, 10 <= z <= 75 km; density weighting for centered field norms','streamfunction_definition':'rho*v=-dPsi/dz; rho*w=dPsi/dy; kg/(m s); Eulerian zonal mean','limitations':['Centered resolved covariance budget, not exact discrete momentum closure. Residual includes viscosity, numerical stresses and diagnostic errors.','No test of residual/Lagrangian circulation or thermodynamic indirectness.','Time-window persistence is a descriptive check, not a statistical significance test.','Stable cap, sponge, free-slip rigid walls and imposed thermal feedback can affect circulation.'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
summary['half_window_period_days']=[250,300]
summary['fifty_day_comparison']=[{'days':windows[i][0],'streamfunction_cosine_vs_previous':cosflat(psi[i-3],psi[i]) if i>2 else None,'zonal_u_cosine_vs_previous':cosine(means[i-3]['zonal_u'],means[i]['zonal_u']) if i>2 else None} for i in [2,5,8,11,14,17]]
summary['reproduction_check']={}
old=json.loads(Path('rendered/ferrel35-day200-audit-v2/summary.json').read_text())
current={r['window_earth_days']:r for r in reports}
checks=0
def check_metrics(a,b,label):
 global checks
 for k,v in a.items():
  if isinstance(v,dict):check_metrics(v,b[k],label+'.'+k)
  elif isinstance(v,(int,float)):
   assert np.isclose(v,b[k],rtol=1e-12,atol=1e-14),(label,k,v,b[k])
   checks+=1
for r in old['windows']:check_metrics(r,current[r['window_earth_days']],r['window_earth_days'])
summary['reproduction_check']={'earlier_day200_all_window_metrics':'PASS','scalar_checks':checks,'rtol':1e-12,'atol':1e-14,'reference_sha256':hashlib.sha256(Path('rendered/ferrel35-day200-audit-v2/summary.json').read_bytes()).hexdigest()}
summary['half_window_comparisons']=[{'days':[start,start+50],'streamfunction_cosine':cosflat(psi[i],psi[i+1]),'smoothed_streamfunction_cosine':cosflat(smooth(psi[i]),smooth(psi[i+1])),'v_cosine':cosine(means[i]['zonal_v'],means[i+1]['zonal_v']),'w_cosine':cosine(means[i]['zonal_w'],means[i+1]['zonal_w'])} for start,i in zip([0,50,100,150,200,250],[0,3,6,9,12,15])]
summary['whole_run_mean_note']='0–300 includes initial adjustment; it is not treated as a stationary circulation estimate.'
summary['numerical_budget_followup']={'source':'rendered/ferrel35-budget200-210-v4/PASS.json','review':'data/ferrel35-width-v4/diagnostic/segment_1/native_budget_review.json','scope':'Separate 200–210 restart only, not an exact budget for the full 300-day run.'}
summary['limitations'].append('The separate 200–210 native-budget restart finds a substantial native-minus-centered advection remainder. The centered residual below is not a measurement of an omitted physical force; its interpretation cannot be repaired by longer averaging alone.')
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
def make_figure(indices,filename,ncols=3):
    groups=(len(indices)+ncols-1)//ncols
    fig,axs=plt.subplots(2*groups,ncols,figsize=(5*ncols,7.5*groups),layout='constrained',squeeze=False)
    ps=[smooth(psi[i]) for i in indices];lim=max(abs(p).max() for p in ps)/1e5
    ul=max(abs(smooth(means[i]['zonal_u'])).max() for i in indices)
    for j,(i,p) in enumerate(zip(indices,ps)):
     c=j%ncols;r=2*(j//ncols);label=windows[i][0];a=means[i]
     im=axs[r,c].pcolormesh((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
     axs[r,c].contour((np.arange(ny+1)*dy-m['ly']/2)/1000,np.arange(nz+1)*dz/1000,p/1e5,levels=np.linspace(-lim,lim,13),colors='k',linewidths=.35)
     axs[r,c].set(title=f'Days {label}: mass streamfunction',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[r,c],label='Psi / 1e5 [kg m^-1 s^-1]')
     im=axs[r+1,c].pcolormesh(y/1000,z/1000,smooth(a['zonal_u']),cmap='RdBu_r',vmin=-ul,vmax=ul,shading='auto');axs[r+1,c].set(title=f'Days {label}: zonal velocity',xlabel='y north [km]',ylabel='Height [km]');fig.colorbar(im,ax=axs[r+1,c],label='u [m/s]')
     for ax in axs[r:r+2,c]:ax.axhline(80,color='k',ls='--',lw=.6)
    fig.suptitle('35N heated beta-plane | time-mean overturning audit\nSmoothing: 200 km in y, 5 km in z; shared limits for each field; dashed line = initial cap base')
    fig.savefig(out/filename,dpi=130);plt.close(fig)
make_figure([15,16,17],'overturning_late.png')
make_figure([2,5,8,11,14,17],'overturning_evolution.png')
make_figure([18],'overturning_whole_run.png',ncols=1)
print(json.dumps({'status':summary['status'],'frames':summary['frames_analyzed'],'reproduction_check':summary['reproduction_check'],'last50':reports[17],'late_persistence':summary['half_window_comparisons'][-1]},indent=2),flush=True)

"""Point-4 hypothesis tests. Run only on Unity, from saved Float64 3D velocities.
Centered continuum reconstructions are not the native WENO momentum budget.
"""
from pathlib import Path
import json,sys,tomllib,hashlib
import numpy as np
import h5py
from scipy.fft import fft,ifft,dst,idst
from scipy.ndimage import gaussian_filter, distance_transform_edt
from scipy.optimize import least_squares
from scipy.special import expit

def dx(q,h):return (np.roll(q,-1,0)-np.roll(q,1,0))/(2*h)
def dy(q,h,parity=1):
 p=np.concatenate((parity*q[:,0:1,...],q,parity*q[:,-1:,...]),axis=1)
 return (p[:,2:,...]-p[:,:-2,...])/(2*h)
def dz(q,h):
 p=np.concatenate((q[:,:,0:1],q,q[:,:,-1:]),axis=2)
 return (p[:,:,2:]-p[:,:,:-2])/(2*h)
def poisson(z,hx,hy):
 # Cell-centered Dirichlet streamfunction at y walls, periodic x.
 nx,ny=z.shape;kx=np.arange(nx);ky=np.arange(1,ny+1)
 lam=-4*np.sin(np.pi*kx/nx)**2/hx**2
 lam=lam[:,None]-4*np.sin(np.pi*ky/(2*ny))[None,:]**2/hy**2
 return idst(ifft(dst(fft(z,axis=0),type=2,axis=1,norm='ortho')/lam,axis=0).real,type=2,axis=1,norm='ortho')
def lap(q,hx,hy):
 p=np.concatenate((-q[:,0:1],q,-q[:,-1:]),axis=1)
 return (np.roll(q,-1,0)-2*q+np.roll(q,1,0))/hx**2+(p[:,2:]-2*q+p[:,:-2])/hy**2

def smooth(q,km,hx,hy):
 # Gaussian FWHM, not a sharp Fourier wavelength cutoff. Odd curl extension at walls.
 if not km:return q.copy()
 sig=km*1000/np.sqrt(8*np.log(2))
 ny=q.shape[1]; ext=np.concatenate((q,-q[:,::-1]),axis=1)
 return gaussian_filter(ext,(sig/hx,sig/hy),mode='wrap')[:,:ny]
def moments(q,weights=None):
 q=np.asarray(q).ravel();w=np.ones_like(q)/len(q) if weights is None else np.broadcast_to(weights,np.shape(weights)).ravel().copy()
 if len(w)!=len(q):raise ValueError('Moment weights have wrong shape')
 w=w/w.sum();m=float(w@q);s2=float(w@(q*q));c2=float(w@((q-m)**2))
 return {'mean':m,'rms':np.sqrt(s2),'raw_third_moment':float(w@(q**3)),
  'raw_skewness':float(w@(q**3)/s2**1.5) if s2>1e-40 else None,
  'central_skewness':float(w@((q-m)**3)/c2**1.5) if c2>1e-40 else None,
  'positive_fraction':float(w@(q>0)),'negative_fraction':float(w@(q<0))}

def field_analysis(p):
 with h5py.File(p) as d:
  a={k:np.asarray(d[k]).T for k in ['u_face','v_face','w_face']}
  rho=d['density'][:];hx=float(d['dx'][()]);hy=float(d['dy'][()]);hz=float(d['dz'][()]);f=float(d['f'][()]);t=float(d['t'][()])
 uf,vf,wf=a['u_face'],a['v_face'],a['w_face']
 u=(uf+np.roll(uf,-1,0))/2;v=(vf[:,:-1,:]+vf[:,1:,:])/2;w=(wf[:,:,:-1]+wf[:,:,1:])/2
 weights=rho[None,None,:]/rho.sum(); avg=lambda q:np.sum(weights*q,axis=2)
 U,Vbar=avg(u),avg(v);zeta=dx(v,hx)-dy(u,hy);Z=avg(zeta)
 up=u-U[:,:,None];vp=v-Vbar[:,:,None];zp=zeta-Z[:,:,None]
 delta=dx(up,hx)+dy(vp,hy,-1)
 sourceV=-2*avg(zp*delta)
 sourceH=-avg(up*dx(zp,hx)+vp*dy(zp,hy,-1))+avg((dy(w,hy)-dz(vp,hz))*dx(w,hx)+(dz(up,hz)-dx(w,hx))*dy(w,hy))
 pred=2*(f+Z)*avg(delta**2)
 adv=-U*dx(Z,hx)-Vbar*dy(Z,hy,-1)
 divU=dx(U,hx)+dy(Vbar,hy,-1);comp=-(2*Z+f)*divU
 nx,ny=Z.shape;y=(np.arange(ny)+.5)*hy-ny*hy/2
 psi=poisson(Z,hx,hy);psi-=np.mean(U)*(y[None,:]+ny*hy/2);psi-=psi.mean()
 Zeddy=Z-Z.mean(0,keepdims=True);psieddy=poisson(Zeddy,hx,hy);psieddy-=psieddy.mean()
 native_div=(np.roll(uf,-1,0)-uf)/hx+np.diff(vf,axis=1)/hy
 # Exact density-weighted horizontal divergence uses native faces; column cancellation must hold.
 divmass=avg(native_div)
 row={'time_days':t/86400,'f':f,'column_mass_kg_m2':float(rho.sum()*hz),'density_kg_m3':rho.tolist(),
      'grid':[nx,ny,len(rho)],'dx_km':hx/1000,'native_depth_mean_divergence_max':float(abs(divmass).max()),
      'centered_depth_mean_divergence_max':float(abs(divU).max()),'moments':{},'source':{}}
 masks={'full':np.ones(Z.shape,dtype=bool),'interior':np.broadcast_to(abs(y)[None,:]<=.4*ny*hy,Z.shape),'wall_strips':np.broadcast_to(abs(y)[None,:]>.4*ny*hy,Z.shape)}
 for region,mask in masks.items():
  row['moments'][region]={}
  for kind,q in [('full_Z',Z),('nonzonal_Z',Zeddy),('zonal_Z',np.broadcast_to(Z.mean(0),Z.shape))]:
   row['moments'][region][kind]={}
   for scale in [0,50,100,200]:
    filtered=smooth(q,scale,hx,hy)
    m=moments(filtered[mask]);hist,bins=np.histogram(filtered[mask]/abs(f),bins=100,density=True)
    m['PDF_zeta_over_abs_f']={'bin_edges':bins.tolist(),'density':hist.tolist()}
    m['filter_FWHM_km']=scale;m['filter_FWHM_cells']=scale*1000/hx
    m['filter_resolved_at_4_cells']=scale==0 or scale*1000/hx>=4
    m['cyclonic_fraction']=float(np.mean(filtered[mask]*f>0))
    row['moments'][region][kind][str(scale)]=m
  zz=zeta[mask,:];ww=np.broadcast_to(weights,zeta.shape)[mask,:]
  volum=moments(zz);density=moments(zz,ww)
  row['moments'][region]['three_dimensional_uniform_volume']=volum
  row['moments'][region]['three_dimensional_density_weighted']=density
  row['moments'][region]['FuSun_reference']={'Ro_g_volume':volum['rms']/abs(f),'2p5_Ro_volume':2.5*volum['rms']/abs(f),
      'rotation_signed_skewness_volume':None if volum['raw_skewness'] is None else np.sign(f)*volum['raw_skewness'],
      'caution':'Fu–Sun relation is a finite-amplitude rotating Boussinesq convection result for 3D vorticity; applicability to this anelastic capped channel is being tested, not assumed.'}
  rms=lambda a:float(np.sqrt(np.mean(a[mask]**2)))
  w2=avg(w*w);active=mask & (w2>=np.quantile(w2[mask],.75)) & (w2>1e-30)
  cores=np.sign(f)*Z>np.sqrt(np.mean(Z[mask]**2))
  # Periodic x distance to the specified cyclonic cores; y is not wrapped.
  tiled=np.tile(cores,(3,1));distance=distance_transform_edt(~tiled,sampling=(hx,hy))[nx:2*nx]
  ring=mask & (~cores) & (distance<=200000)
  conditional=lambda a,mm:float(a[mm].mean()) if mm.any() else None
  x=pred[mask];yv=sourceV[mask];den=float(x@x)
  tau=float(x@yv/den) if den>1e-70 else None
  corr=float(np.corrcoef(x,yv)[0,1]) if x.std()>0 and yv.std()>0 else None
  ss={'V_mean':float(sourceV[mask].mean()),'H_mean':float(sourceH[mask].mean()),'V_plus_H_mean':float((sourceV+sourceH)[mask].mean()),
      'V_in_active_convection_mean':conditional(sourceV,active),'rotation_signed_V_in_active_convection_mean':conditional(np.sign(f)*sourceV,active),'H_in_200km_core_ring_mean':conditional(sourceH,ring),'rotation_signed_H_in_200km_core_ring_mean':conditional(np.sign(f)*sourceH,ring),'core_definition':'sign(f) Z > interior/domain RMS Z; ring is outside these cores within 200 km','V_rms':rms(sourceV),'H_rms':rms(sourceH),'V_plus_H_rms':rms(sourceV+sourceH),'V_positive_fraction':float(np.mean(sourceV[mask]>0)),
      'effective_tau_seconds_through_origin':tau,'tau_predictor_correlation':corr,
      'tau_fit_R2':None if tau is None or yv.var()==0 else float(1-np.mean((yv-tau*x)**2)/yv.var()),'third_moment_source_s_inv4':{}}
  for name,q in [('V',sourceV),('H',sourceH),('barotropic_advection',adv),('mean_divergence',comp)]:
   ss['third_moment_source_s_inv4'][name]=float(np.mean(3*Z[mask]**2*q[mask]))
   low=smooth(q,200,hx,hy);high=q-low
   ss[name+'_200km_split']={'low_rms':rms(low),'high_rms':rms(high),'low_third_moment_source':float(np.mean(3*Z[mask]**2*low[mask])),'high_third_moment_source':float(np.mean(3*Z[mask]**2*high[mask]))}
  row['source'][region]=ss
 # Reduce volume only after all 3D derivatives and density-weighted products have been evaluated.
 maps={'Z':Z,'Psi':psi,'Z_nonzonal':Zeddy,'Psi_nonzonal':psieddy,'U':U,'Vbar':Vbar,'vertical_source':sourceV,'horizontal_source':sourceH,
       'barotropic_advection':adv,'mean_divergence_source':comp,'tau_predictor':pred,'convective_w2':avg(w*w),'x':(np.arange(nx)+.5)*hx,'y':y,'time_days':t/86400}
 return row,maps

def fit_relation(x,y,seed=934):
 x=np.asarray(x);y=np.asarray(y);xc=x.mean();xs=x.std();yc=y.mean();ys=y.std()
 if xs<=1e-25 or ys<=1e-25:return {'status':'DEGENERATE'}
 q=(x-xc)/xs;r=(y-yc)/ys
 # Equal-spatial thinning only for optimization cost; evaluate R² on all retained pooled samples.
 index=np.linspace(0,len(q)-1,min(3000,len(q))).astype(int);qt=q[index];rt=r[index]
 co=np.polynomial.polynomial.polyfit(qt,rt,2);quad=np.polynomial.polynomial.polyval(q,co)
 raw_c=ys*co[2]/xs**2;raw_b=ys*(co[1]/xs-2*co[2]*xc/xs**2)
 def sigmoid(par,q):return par[0]+np.exp(par[1])*expit(par[2]+par[3]*q)
 starts=[]
 for slope in [-2.,2.]:
  fit=least_squares(lambda a:sigmoid(a,qt)-rt,[-2.,np.log(4.),0.,slope],bounds=([-30,-8,-30,-30],[30,8,30,30]),max_nfev=150)
  starts.append(fit)
 fit=min(starts,key=lambda a:np.mean(a.fun**2));par=fit.x;prediction=sigmoid(par,q)
 sigma2=yc+ys*par[0];sigma1=sigma2+ys*np.exp(par[1]);occupancy=(yc-sigma2)/(sigma1-sigma2)
 r2=lambda pred:float(1-np.mean((r-pred)**2)/np.var(r))
 return {'status':'FIT','samples':len(q),'Psi_mean':xc,'Psi_standard_deviation':xs,'Z_mean':yc,'Z_standard_deviation':ys,
   'quadratic_standardized_coefficients':co.tolist(),'quadratic_R2':r2(quad),'raw_b':raw_b,'raw_c':raw_c,
   'raw_c_over_b':None if abs(raw_b)<1e-10*abs(ys/xs) else raw_c/raw_b,
   'standardized_c_over_b':None if abs(co[1])<1e-10 else float(co[2]/co[1]),
   'logistic_R2':r2(prediction),'logistic_sigma1':sigma1,'logistic_sigma2':sigma2,
   'logistic_standardized_intercept':float(par[2]),'logistic_standardized_slope':float(par[3]),
   'logistic_inferred_occupancy_from_mean':occupancy,'logistic_parameters_at_bounds':bool(np.any(abs(par-np.array([-30,-8,-30,-30]))<1e-3) or np.any(abs(par-np.array([30,8,30,30]))<1e-3)),
   'identifiability':'a and mu cannot be separately identified from this curve; only a combined intercept is fitted. Occupancy inferred from the mean is not an independently fitted cyclone area.',
   'equilibrium_warning':'Explained variance is association, not proof of stationarity or third-invariant conservation.'}

def analyze(run,out):
 out.mkdir(parents=True,exist_ok=False);meta=tomllib.loads((run/'metadata.toml').read_text())
 assert meta['status']=='COMPLETE'
 rows=[];samples={k:[] for k in ['full','interior','nonzonal_interior']};previous=None;budgets=[]
 files=sorted((run/'fields3d').glob('state_*.jld2'));assert len(files)>=2
 for p in files:
  row,m=field_analysis(p);rows.append(row);t=m['time_days'];mask=np.broadcast_to(abs(m['y'])[None,:]<=.4*meta['ly'],m['Z'].shape)
  for name,zz,pp,mm in [('full',m['Z'],m['Psi'],np.ones_like(mask)),('interior',m['Z'],m['Psi'],mask),('nonzonal_interior',m['Z_nonzonal'],m['Psi_nonzonal'],mask)]:
   # Keep time-balanced deterministic spatial samples for plotting and pooled fits.
   ix=np.flatnonzero(mm.ravel())[::max(1,int(mm.sum()/12000))]
   samples[name].append(np.column_stack((pp.ravel()[ix],zz.ravel()[ix],np.full(len(ix),t))))
   row.setdefault('fits',{})[name]=fit_relation(pp[mm],zz[mm])
   zzm=zz[mm];rms=np.sqrt(np.mean(zzm**2))
   row['fits'][name]['measured_core_fraction_zeta_over_rms_gt1']=float(np.mean(np.sign(row['f'])*zzm>rms))
   row['fits'][name]['measured_positive_core_fraction_zeta_over_rms_gt1']=float(np.mean(zzm>rms))
   row['fits'][name]['measured_cyclonic_fraction']=float(np.mean(zzm*row['f']>0))
  np.savez_compressed(out/(p.stem+'_maps.npz'),**m)
  if previous is not None:
   pr,pm=previous;dt=(t-pm['time_days'])*86400
   storage=(m['Z']-pm['Z'])/dt
   rhs=sum((m[k]+pm[k])/2 for k in ['vertical_source','horizontal_source','barotropic_advection','mean_divergence_source'])
   remainder=storage-rhs
   budgets.append({'start_day':pm['time_days'],'end_day':t,'storage_rms':float(np.sqrt(np.mean(storage**2))),
       'centered_inviscid_reconstruction_rms':float(np.sqrt(np.mean(rhs**2))),'unclosed_remainder_rms':float(np.sqrt(np.mean(remainder**2))),
       'third_moment_storage':float(np.mean(m['Z']**3-pm['Z']**3)/dt),
       'third_moment_reconstructed_tendency':float(sum((row['source']['full']['third_moment_source_s_inv4'][k]+pr['source']['full']['third_moment_source_s_inv4'][k])/2 for k in ['V','H','barotropic_advection','mean_divergence'])),
       'meaning':'Remainder includes viscosity, discrete reconstruction, temporal sampling and omitted terms; this is not a closed native budget.'})
  previous=(row,m)
  print('ANALYZED',p,t,flush=True)
 pooled={}
 for name,ss in samples.items():
  a=np.concatenate(ss);np.savez_compressed(out/(name+'_scatter.npz'),Psi=a[:,0],Z=a[:,1],day=a[:,2])
  pooled[name]=fit_relation(a[:,0],a[:,1])
 report={'status':'COMPLETE','compute_location':'Unity','source_run':str(run),'metadata':meta,'snapshot_count':len(files),'snapshots':rows,'pooled_fits':pooled,'finite_interval_budgets':budgets,
   'limitations':['Dry anelastic channel, not a validation of a moist or spherical Jupiter model.','50 km FWHM spans only 2.13 coarse cells and is unresolved under the four-cell screening criterion.',
    'Two-level occupancy and offset are degenerate; no independent free a and mu fit is reported.','The Fu–Sun scaling concerns 3D vorticity; barotropic and 3D statistics are kept separate.',
    'V and H are continuum centered reconstructions; budget remainders are retained rather than called dissipation.',
    'A five-day from-rest run may not develop a condensate; absence at day 5 is not evidence against eventual sign selection.',
    'Flipping f changes the sign of cyclone vorticity, not its cyclonic classification.',
    'All reference-density values and column masses are recorded; masks separate the outer 10 percent at each channel wall and zonal mean.'],
   'analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 def clean(x):
  if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
  if isinstance(x,(list,tuple)):return [clean(v) for v in x]
  if isinstance(x,np.generic):return clean(x.item())
  if isinstance(x,float) and not np.isfinite(x):return None
  return x
 (out/'diagnostics.json').write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n')
 return report

if __name__=='__main__':analyze(Path(sys.argv[1]),Path(sys.argv[2]))

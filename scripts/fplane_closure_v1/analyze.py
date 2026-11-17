"""Existing-field closure test. Execute only on Unity; never advances a model."""
from pathlib import Path
import sys, json, csv, tomllib, hashlib, time
import numpy as np
import h5py
from scipy.ndimage import gaussian_filter
from scipy.optimize import lsq_linear
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FEATURES=['P0(z)', 'P0(h2_laplacian_z)', 'sign_f_P0(z2)', '-P0(z3)',
          'sign_f_P0(strain2)', 'P0(z_strain2)', 'sign_f_P0(h2_grad_z_squared)',
          'P0(z_convective_intensity)', 'sign_f_P0(h2_laplacian_z2)']
MODELS={'linear':[0,1], 'polynomial_free':[0,1,2,3],
        'polynomial_sign_saturating':[0,1,2,3], 'augmented':list(range(9))}
TERMS=['mean_transport','mean_flux_minus_transport','convective_stress',
       'native_advection_minus_centered_flux','native_coriolis','native_diffusion','native_horizontal_forcing']

def dx(q,h): return (np.roll(q,-1,axis=0)-np.roll(q,1,axis=0))/(2*h)
def dy(q,h,parity=1):
    p=np.concatenate((parity*q[:,0:1],q,parity*q[:,-1:]),axis=1)
    return (p[:,2:]-p[:,:-2])/(2*h)
def curl(u,v,h): return dx(v,h)-dy(u,h)
def stress(xx,xy,yy,h): return curl(-dx(xx,h)-dy(xy,h,-1),-dx(xy,h)-dy(yy,h),h)
def lap(q,h,parity=-1): return dx(dx(q,h),h)+dy(dy(q,h,parity),h,-parity)
def p0(q): return q-q.mean()
def smooth(q,sigma,parity):
    if sigma==0:return q
    n=q.shape[1]
    ext=np.concatenate((q,parity*q[:,::-1]),axis=1)
    return gaussian_filter(ext,sigma,mode='wrap')[:,:n]
def read(p,keys):
    with h5py.File(p) as h:return {k:np.asarray(h[k]).T for k in keys}
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.bool_):return bool(x)
    return x
def save(p,x):p.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n')

def features(U,V,Rxx,Rxy,Ryy,f,h,sigma=0):
    # Filter the convective stress itself. Barotropic subfilter stresses are NOT folded into C.
    U=smooth(U,sigma,1);V=smooth(V,sigma,-1)
    xx=smooth(Rxx,sigma,1);xy=smooth(Rxy,sigma,-1);yy=smooth(Ryy,sigma,1)
    Z=curl(U,V,h);C=stress(xx,xy,yy,h);F=abs(f);z=Z/F
    ux,uy,vx,vy=dx(U,h)/F,dy(U,h)/F,dx(V,h)/F,dy(V,h,-1)/F
    s2=(ux-vy)**2+(uy+vx)**2
    intensity=(xx+yy)/(F*1.2e6)**2
    grad2=h*h*(dx(z,h)**2+dy(z,h,-1)**2)
    raw=[z,h*h*lap(z,h),np.sign(f)*z*z,-z**3,np.sign(f)*s2,z*s2,
         np.sign(f)*grad2,z*intensity,np.sign(f)*h*h*lap(z*z,h,1)]
    X=np.stack([p0(q) for q in raw],axis=-1)
    assert np.isfinite(X).all() and np.isfinite(C).all()
    return X,C/(F*F),z,U,V

def low_modes(n,h):
    # Orthogonal channel streamfunctions, odd at y walls. Matched centered derivative symbols.
    L=n*h;x=(np.arange(n)+.5)*h;y=x.copy();modes=[];lam=[]
    for m in range(3):
        for k in range(1,5):
            if m*m+(k/2)**2>4+1e-12:continue
            sx=np.sin(2*np.pi*m/n)/h;sy=np.sin(np.pi*k/n)/h
            for trig in ([np.cos] if m==0 else [np.cos,np.sin]):
                q=trig(2*np.pi*m*x[:,None]/L)*np.sin(np.pi*k*y[None,:]/L)
                q/=np.sqrt(np.mean(q*q));modes.append(q);lam.append(sx*sx+sy*sy)
    return np.stack(modes),np.asarray(lam)

def fit_model(G,r,indices,constrained=False):
    G=G[np.ix_(indices,indices)];r=r[indices]
    scale=np.sqrt(np.maximum(np.diag(G),1e-300));A=G/scale[:,None]/scale[None,:]
    eig,V=np.linalg.eigh(A);keep=eig>max(eig.max(),1)*1e-11
    assert eig.min()>-1e-9 and keep.sum()>0
    design=np.sqrt(eig[keep])[:,None]*V[:,keep].T
    target=(V[:,keep].T@(r/scale))/np.sqrt(eig[keep])
    if constrained:
        lower=np.array([0. if j in [2,3] else -np.inf for j in indices])
        sol=lsq_linear(design,target,bounds=(lower,np.inf),tol=1e-12,max_iter=500)
        assert sol.success,sol.message
        c=sol.x/scale
    else:c=np.linalg.lstsq(design,target,rcond=1e-11)[0]/scale
    out=np.zeros(9);out[indices]=c
    return out,{'rank':int(keep.sum()),'columns':len(indices),'scaled_condition':float(eig[keep].max()/eig[keep].min())}

def validation():
    n=64;h=2*np.pi/n;x=(np.arange(n)+.5)*h;y=x.copy()
    psi=np.cos(2*x[:,None])*np.sin(y[None,:]/2)
    U=-dy(psi,h,-1);V=dx(psi,h);Z=curl(U,V,h)
    lam=(np.sin(2*h)/h)**2+(np.sin(h/2)/h)**2
    err=np.max(abs(Z+lam*psi));assert err<1e-11
    Rxx=2+.2*np.cos(x[:,None])*np.cos(y[None,:]/2)
    Ryy=1+.3*np.sin(2*x[:,None])*np.cos(y[None,:])
    Rxy=.1*np.cos(x[:,None])*np.sin(y[None,:]/2)
    for sigma in [0,4]:
        X,C,z,*_=features(U,V,Rxx,Rxy,Ryy,1,h,sigma)
        Xm,Cm,zm,*_=features(-U[::-1],V[::-1],Rxx[::-1],-Rxy[::-1],Ryy[::-1],-1,h,sigma)
        assert np.max(abs(Cm+C[::-1]))<1e-11
        assert np.max(abs(Xm+X[::-1]))<1e-9
        assert np.max(abs(zm+z[::-1]))<1e-11
        assert np.max(abs(X.mean((0,1))))<1e-11
    rng=np.random.default_rng(142);D=rng.normal(size=(1000,9));coef=np.array([.1,-.2,.3,.4,0,0,0,0,0]);target=D@coef
    for constraint in [False,True]:
        got,_=fit_model(D.T@D,D.T@target,[0,1,2,3],constraint)
        assert np.max(abs(got-coef))<1e-10
    z=rng.normal(size=(64,64));b=.2
    actual=3*np.mean(z*z*b*p0(z*z));expected=3*b*(np.mean(z**4)-np.mean(z*z)**2)
    assert abs(actual-expected)<1e-12
    modes,lamb=low_modes(n,h)
    ortho=modes.reshape(len(modes),-1)@modes.reshape(len(modes),-1).T/(n*n)
    assert np.max(abs(ortho-np.eye(len(modes))))<1e-12
    # Work from curl equals work from the velocity force for a resolved mode.
    phi=modes[0];fu=np.sin(x[:,None])*np.cos(y[None,:]/2)+.3*U;fv=.2*V
    lhs=np.mean((-dy(phi,h,-1))*fu+dx(phi,h)*fv)
    rhs=-np.mean(phi*curl(fu,fv,h));assert abs(lhs-rhs)<1e-12
    return {'status':'PASS','analytic_centered_curl_error':err,'reflection_all_features_both_filters':True,
            'synthetic_coefficient_recovery':True,'quadratic_third_moment_identity':True,
            'low_mode_orthogonality_and_work_identity':True,'no_model_evolution':True}

def record(path,meta,sigma):
    d=read(path,['t','ubar','vbar','Rxx','Rxy','Ryy'])
    f=float(meta['f']);h=meta['lx']/meta['nx'];assert meta['lx']/meta['nx']==meta['ly']/meta['ny']
    X,q,z,U,V=features(*(d[k] for k in ['ubar','vbar','Rxx','Rxy','Ryy']),f,h,sigma)
    return {'X':X,'q':q,'z':z,'U':U,'V':V,'t':float(d['t']),'f':f,'h':h}

def scores(y,p):
    e=p-y;den=np.mean(y*y);var=np.var(y)
    return {'nrmse':np.sqrt(np.mean(e*e)/max(den,1e-300)),
            'r2':1-np.mean(e*e)/max(var,1e-300),
            'correlation':np.corrcoef(y,p)[0,1] if np.std(p)>0 and np.std(y)>0 else np.nan,
            'observed_mean':np.mean(y),'predicted_mean':np.mean(p),'source_rms':np.sqrt(den)}

def main():
    root,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False);started=time.time()
    valid=validation();save(out/'validation.json',valid)
    src=root/'fplane-point4-v2/late_heated';meta=tomllib.loads((src/'metadata.toml').read_text())
    assert meta['status']=='COMPLETE' and meta['rotation_model']=='f-plane' and meta['beta']==0
    assert meta['heat_flux']==162 and [meta[k] for k in ['nx','ny','nz']]==[256,256,64]
    n=meta['nx'];h=meta['lx']/n;f=abs(meta['f']);modes,lamb=low_modes(n,h)
    ym=abs((np.arange(n)+.5)/n-.5)<=.4;mask=np.broadcast_to(ym[None,:],(n,n));wall=~mask
    files=sorted(src.glob('frame_*.jld2'));assert len(files)==81
    report={'status':'COMPLETE','compute_location':'Unity','new_simulations':False,
      'equation':'dZ/dt + U.grad Z = C + mean-centering correction + native-minus-centered advection + Coriolis + diffusion + horizontal forcing',
      'target':'C=-curl(div R), R_ij=<u_i prime u_j prime>_rho; centered diagnostic kept distinct from native WENO',
      'physical_setup':{k:meta[k] for k in ['nx','ny','nz','lx','ly','height','f','beta','heat_flux','viscosity','cap_height']},
      'normalization':'z=Z/abs(f), q=C/abs(f)^2; tprime=abs(f)*t. Feature derivatives use h=grid spacing. Convective intensity=(Rxx+Ryy)/(abs(f)*1200km)^2.',
      'features':FEATURES,'training_days':[90,95],'training_interval':'90 <= day < 95',
      'test_days':[95,100],'test_interval':'95 <= day <= 100','fit_region':'central 80% of north-south channel',
      'filters':{},'limits':['No new simulation or changed cap/forcing.','Fit coefficients are inferred, not copied from dimensionless 2D tests.',
      'Spatial samples are correlated; no independent-sample significance claim.','Time-block separation does not supply an independent realization.',
      'All feature means are removed over the full channel. Actual target means are retained and reported; channel circulation is not assumed zero.',
      'Extra linear z and Laplacian terms allow a more flexible linear response than the arbitrary 2D forcing band; this does not derive that band.',
      'Gaussian filtering is a diagnostic sensitivity, not grid refinement; it filters the same convective stress and does not add barotropic subfilter stresses.',
      'Snapshot source work is not a closed ten-day energy or third-moment budget. Native six-hour accounting is reported separately.',
      'The f-minus audit inherits the +f mature state and contains adjustment; it is not a mirrored independent realization.']}
    all_rows=[];predictions={};conditional={};plot_cache={};coefs={}
    for sigma in [0,4]:
        tag=f'sigma{sigma}';records=[record(p,meta,sigma) for p in files]
        days=np.array([r['t']/86400 for r in records]);assert np.allclose(days,np.linspace(90,100,81))
        train=days<95;test=~train
        G=np.zeros((9,9));rhs=np.zeros(9)
        for r,use in zip(records,train):
            if use:
                X=r['X'][mask];y=r['q'][mask];G+=X.T@X;rhs+=X.T@y
        coeff={};details={}
        for name,idx in MODELS.items():coeff[name],details[name]=fit_model(G,rhs,idx,name=='polynomial_sign_saturating')
        coefs[tag]=coeff
        y=np.concatenate([r['q'][mask] for r,use in zip(records,test) if use])
        models={};preds={}
        for name,c in coeff.items():
            pp=[r['X']@c for r in records];preds[name]=pp
            p=np.concatenate([a[mask] for a,use in zip(pp,test) if use])
            models[name]={'coefficients':dict(zip(FEATURES,c)),'fit':details[name],'heldout_interior':scores(y,p)}
            for region,sel in [('full',np.ones_like(mask)),('walls',wall)]:
                models[name]['heldout_'+region]=scores(np.concatenate([r['q'][sel] for r,use in zip(records,test) if use]),np.concatenate([a[sel] for a,use in zip(pp,test) if use]))
        for it,r in enumerate(records):
            z=r['z'];Z=z*f;zc=z[mask]-np.mean(z[mask]);ql=r['q'];zm=np.mean(modes*Z,axis=(1,2));psi=-np.sum((zm/lamb)[:,None,None]*modes,axis=0)
            # Independent modal check from velocity projections, avoiding an arbitrary Poisson inversion.
            for mode,lam,zz in zip(modes,lamb,zm):
                velocity_projection=np.mean(-dy(mode,h,-1)*r['U']+dx(mode,h)*r['V'])
                assert abs(velocity_projection+zz)<1e-12*max(1,abs(zz))
            for name,p in [('observed',ql)]+[(name,arr[it]) for name,arr in preds.items()]:
                row={'sigma_cells':sigma,'day':days[it],'split':'train' if train[it] else 'test','model':name,
                  'central_m3_source_dimless':3*np.mean(zc*zc*(p[mask]-np.mean(p[mask]))),
                  'raw_m3_source_dimless':3*np.mean(z[mask]**2*p[mask]),
                  'low_mode_energy_work_m2_s3':-np.mean(psi*p*f*f),'low_mode_energy_m2_s2':.5*np.sum(zm*zm/lamb),
                  'source_mean_s_inv2':np.mean(p)*f*f}
                if name!='observed':row.update(scores(ql[mask],p[mask]))
                all_rows.append(row)
        edges=np.linspace(-.5,.5,81);counts=np.zeros(80);zzsum=np.zeros(80);sums={name:np.zeros(80) for name in ['observed']+list(coeff)}
        for it,r in enumerate(records):
            if not test[it]:continue
            z=r['z'][mask];ii=np.digitize(z,edges)-1;ok=(ii>=0)&(ii<80)
            counts+=np.bincount(ii[ok],minlength=80);zzsum+=np.bincount(ii[ok],weights=z[ok],minlength=80)
            for name,a in [('observed',r['q'])]+[(name,arr[it]) for name,arr in preds.items()]:sums[name]+=np.bincount(ii[ok],weights=a[mask][ok],minlength=80)
        conditional[tag]={'edges':edges,'count':counts,'mean_z':np.divide(zzsum,counts,out=np.full(80,np.nan),where=counts>0),
           'means':{name:np.divide(s,counts,out=np.full(80,np.nan),where=counts>0) for name,s in sums.items()},
           'outside_plot_range_fraction':1-counts.sum()/(test.sum()*mask.sum())}
        # Primary science comparisons are temporal held-out errors, with no fit using held-out targets.
        report['filters'][tag]={'sigma_cells':sigma,'sigma_km':sigma*h/1000,'training_frames':int(train.sum()),'test_frames':int(test.sum()),'models':models,
          'input_ranges':{'train_z':[min(r['z'][mask].min() for r,u in zip(records,train) if u),max(r['z'][mask].max() for r,u in zip(records,train) if u)],
                          'test_z':[min(r['z'][mask].min() for r,u in zip(records,test) if u),max(r['z'][mask].max() for r,u in zip(records,test) if u)]}}
        plot_cache[tag]={'Z':records[-1]['z'],'C':records[-1]['q'],'polynomial':preds['polynomial_sign_saturating'][-1],'augmented':preds['augmented'][-1]}
        print('FIT_COMPLETE',tag,{k:round(v['heldout_interior']['r2'],4) for k,v in models.items()},flush=True)
        del records,preds

    # Independently verify depth-weighting and covariance against two actual 3D files.
    fullchecks=[]
    for p in [sorted((src/'fields3d').glob('state_*.jld2'))[0],sorted((src/'fields3d').glob('state_*.jld2'))[-1]]:
        d=read(p,['u_face','v_face','density','t']);u=(d['u_face']+np.roll(d['u_face'],-1,axis=0))/2;v=(d['v_face'][:,:-1]+d['v_face'][:,1:])/2
        wt=d['density']/d['density'].sum();U=np.sum(u*wt,axis=2);V=np.sum(v*wt,axis=2)
        a=read(src/f'frame_{round(float(d["t"])/meta["output_seconds"]):05d}.jld2',['ubar','vbar','Rxx','Rxy','Ryy'])
        computed=[U,V,np.sum((u-U[:,:,None])**2*wt,axis=2),np.sum((u-U[:,:,None])*(v-V[:,:,None])*wt,axis=2),np.sum((v-V[:,:,None])**2*wt,axis=2)]
        errors={}
        for k,q in zip(['ubar','vbar','Rxx','Rxy','Ryy'],computed):
            errors[k]=np.max(abs(q-a[k]));assert errors[k]<1e-10*max(np.max(abs(q)),1)
        fullchecks.append({'day':float(d['t'])/86400,'max_errors':errors})
    valid['independent_3d_density_weighted_checks']=fullchecks

    # Reuse the exact every-stage budget. Never relabel the reconstruction difference as convection.
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'fplane_sign_v1'))
    from gate import check_case,read as budget_read
    audits={}
    for case in ['plus','minus']:
        d=root/'fplane-sign-v1'/case;gate=check_case(d);mt=tomllib.loads((d/'metadata.toml').read_text())
        frames=sorted(d.glob('frame_*.jld2'));rr=[record(p,mt,0) for p in frames]
        bb={round(float(a['t'])):a for a in (budget_read(p) for p in sorted((d/'native_vorticity').glob('budget_*.jld2')))}
        sampled=[];native=[];correction=[];pp={name:[] for name in MODELS};closure_max=0
        for a,b in zip(rr[:-1],rr[1:]):
            dt=b['t']-a['t'];before=bb[round(a['t'])];after=bb[round(b['t'])]
            increments=after['cumulative_parts']-before['cumulative_parts'];parts=increments/dt
            closure_max=max(closure_max,np.max(abs((after['Z']-before['Z'])-increments.sum(2))))
            assert np.max(abs(a['z']*abs(mt['f'])-before['Z']))<1e-16
            sampled.append(.5*(a['q']+b['q'])[mask]);native.append(parts[:,:,2][mask]/abs(mt['f'])**2)
            correction.append(parts[:,:,3][mask]/abs(mt['f'])**2)
            X=.5*(a['X']+b['X'])
            for name,c in coefs['sigma0'].items():pp[name].append((X@c)[mask])
        yy=np.concatenate(native);ss=np.concatenate(sampled);cc=np.concatenate(correction)
        first=bb[round(rr[0]['t'])];last=bb[round(rr[-1]['t'])]
        exactparts=last['central_third_parts']-first['central_third_parts']
        conv=exactparts[1,2];unassigned=abs(exactparts[1,1])+abs(exactparts[1,3])
        audits[case]={'native_gate':gate,'interval_increment_closure_max':closure_max,
          'trapezoidal_snapshot_C_vs_stage_integral':scores(yy,ss),
          'native_advection_reconstruction_rms_over_C':np.sqrt(np.mean(cc*cc)/max(np.mean(yy*yy),1e-300)),
          'combined_effective_source_note':'C plus native-minus-centered advection is a numerical effective source, not pure physical convection.',
          'interior_exact_central_third_term_changes':dict(zip(TERMS,exactparts[1])),
          'interior_unassigned_over_convective_moment':unassigned/max(abs(conv),1e-300),
          'model_transfer_against_stage_integrated_C':{name:scores(yy,np.concatenate(p)) for name,p in pp.items()},
          'model_transfer_against_C_plus_numerical_remainder':{name:scores(yy+cc,np.concatenate(p)) for name,p in pp.items()},
          'interpretation':'Positive f overlaps the training trajectory near day90; negative f is inherited-state adjustment, not a mirror test.'}
        print('NATIVE_AUDIT',case,'sampling_error',audits[case]['trapezoidal_snapshot_C_vs_stage_integral']['nrmse'],flush=True)
    report['native_six_hour_audit']=audits
    for tag,entry in report['filters'].items():
        sigma=entry['sigma_cells']
        observed=[r for r in all_rows if r['sigma_cells']==sigma and r['split']=='test' and r['model']=='observed']
        for name,model in entry['models'].items():
            predicted=[r for r in all_rows if r['sigma_cells']==sigma and r['split']=='test' and r['model']==name]
            assert [r['day'] for r in observed]==[r['day'] for r in predicted]
            model['heldout_source_contributions']={}
            for key in ['central_m3_source_dimless','raw_m3_source_dimless','low_mode_energy_work_m2_s3']:
                a=np.array([r[key] for r in observed]);b=np.array([r[key] for r in predicted])
                model['heldout_source_contributions'][key]=dict(scores(a,b),observed_time_mean=a.mean(),predicted_time_mean=b.mean(),same_sign_fraction=np.mean(np.sign(a)==np.sign(b)))
    report['elapsed_seconds']=time.time()-started
    save(out/'report.json',report);save(out/'validation.json',valid);save(out/'conditional.json',conditional)
    cols=list(dict.fromkeys(k for r in all_rows for k in r))
    with (out/'timeseries.csv').open('w') as fp:
        wr=csv.DictWriter(fp,cols);wr.writeheader();wr.writerows(all_rows)
    plot(report,all_rows,conditional,plot_cache,out)
    manifest={str(p):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
    save(out/'private_provenance.json',{'inputs':manifest,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    save(out/'PASS.json',{'status':'COMPLETE','validation':'PASS','assets':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in out.iterdir() if p.name!='private_provenance.json'}})
    print('CLOSURE_ANALYSIS_COMPLETE',time.time()-started,flush=True)

def plot(report,rows,conditional,maps,out):
    colors={'observed':'black','linear':'#808080','polynomial_free':'#c18424','polynomial_sign_saturating':'#176caa','augmented':'#a1357c'}
    labels={'observed':'Measured convective stress','linear':'Linear','polynomial_free':'Polynomial: free signs','polynomial_sign_saturating':'Polynomial: sign + saturation','augmented':'Add strain / gradients / intensity'}
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for j,tag in enumerate(['sigma0','sigma4']):
        cond=conditional[tag];keep=cond['count']>=50
        for name,a in cond['means'].items():axs[0,j].plot(cond['mean_z'][keep],a[keep],label=labels[name],color=colors[name])
        axs[0,j].set(title=f'Held-out conditional drift: {tag}',xlabel='Depth-mean Z / |f|',ylabel='C / |f|²');axs[0,j].axhline(0,color='.6',lw=.5)
        names=list(MODELS);values=[report['filters'][tag]['models'][k]['heldout_interior']['r2'] for k in names]
        axs[1,j].bar(range(4),values,color=[colors[k] for k in names]);axs[1,j].set(xticks=range(4),xticklabels=[labels[k] for k in names],ylabel='Held-out R² (negative = worse than mean)',title='Fit days 90–95; test days 95–100');axs[1,j].tick_params(axis='x',rotation=22);axs[1,j].axhline(0,color='black',lw=.6)
    axs[0,0].legend(fontsize=8);fig.suptitle('Can local polynomial feedback predict the convective vorticity source?')
    fig.savefig(out/'closure_prediction.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for j,sigma in enumerate([0,4]):
        for name in ['observed','linear','polynomial_sign_saturating','augmented']:
            a=[r for r in rows if r['sigma_cells']==sigma and r['model']==name]
            for i,col in enumerate(['central_m3_source_dimless','low_mode_energy_work_m2_s3']):axs[i,j].plot([r['day'] for r in a],[r[col] for r in a],color=colors[name],label=labels[name],lw=1)
        for i in range(2):axs[i,j].axvline(95,color='.5',ls='--');axs[i,j].axhline(0,color='.5',lw=.5);axs[i,j].set(xlabel='Earth day',title=f'Filter sigma = {sigma} cells')
    axs[0,0].set_ylabel('Central-third-moment source (dimensionless)');axs[1,0].set_ylabel('Source work in large channel modes [m²/s³]');axs[0,0].legend(fontsize=8)
    fig.suptitle('Source contributions: predictions versus measured stress (not total budgets)')
    fig.savefig(out/'closure_moments_energy.png',dpi=150);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(14,4.8),layout='constrained')
    for ax,case in zip(axs[:2],['plus','minus']):
        a=report['native_six_hour_audit'][case];parts=a['interior_exact_central_third_term_changes'];ax.bar(np.arange(7),list(parts.values()));ax.set(xticks=np.arange(7),xticklabels=['Transport','Mean correction','Convective','Native correction','Coriolis','Diffusion','Forcing'],title=f'{case}: exact 6-hour central-moment contributions',ylabel='Change [s⁻³]');ax.tick_params(axis='x',rotation=55)
    a=report['native_six_hour_audit'];x=np.arange(2)
    axs[2].bar(x-.17,[a[k]['trapezoidal_snapshot_C_vs_stage_integral']['nrmse'] for k in ['plus','minus']],width=.34,label='Temporal sampling / C')
    axs[2].bar(x+.17,[a[k]['native_advection_reconstruction_rms_over_C'] for k in ['plus','minus']],width=.34,label='Native correction / C')
    axs[2].set(xticks=x,xticklabels=['+f','−f'],ylabel='RMS ratio',title='Sampling and numerical attribution are separate');axs[2].legend(fontsize=8)
    fig.savefig(out/'closure_native_audit.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()

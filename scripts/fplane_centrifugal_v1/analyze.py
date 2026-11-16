"""Postprocess the existing mature heated f-plane run; never evolve a model."""
from pathlib import Path
import csv,hashlib,json,sys,tomllib,time
import numpy as np
import h5py
from scipy.optimize import linear_sum_assignment
from operators import dx,dy,centered,curl,corner_curl,candidates,radial
from check import checks

def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.bool_):return bool(x)
    return x
def save(path,obj):path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')

src,out=map(Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=False)
started=time.time();validation=checks();save(out/'validation.json',validation)
meta=tomllib.loads((src/'metadata.toml').read_text())
assert meta['status']=='COMPLETE' and meta['rotation_model']=='f-plane' and meta.get('beta',0)==0
assert meta['heat_flux']==162 and meta['nx']==256 and meta['ny']==256 and meta['nz']==64
files=sorted((src/'fields3d').glob('state_*.jld2'));assert len(files)==11
frames=sorted(src.glob('frame_*.jld2'));assert len(frames)==81
bins=np.linspace(-4,4,401);mom=[];face=[];hists=[];energy=[];allcores=[];profiles=[];manifest=[];extremes=[]
for p in files:
    with h5py.File(p) as h:
        uf,vf,wf,theta=[np.asarray(h[k]).T for k in ['u_face','v_face','w_face','theta']]
        rho=np.asarray(h['density']);z=np.asarray(h['z']);hx,hy,hz,f,t=[float(h[k][()]) for k in ['dx','dy','dz','f','t']]
    assert f>0 and np.all(np.isfinite(uf)) and np.all(np.isfinite(theta))
    nx,ny,nz=theta.shape;Lx=nx*hx;Ly=ny*hy;y=(np.arange(ny)+.5)*hy-Ly/2
    mask=np.abs(y)<=.4*Ly;wall=~mask
    u,v,w=centered(uf,vf,wf);q=curl(u,v,hx,hy);a=q/f
    native=corner_curl(uf,vf,hx,hy)/f;yf=np.arange(ny+1)*hy-Ly/2;maskf=abs(yf)<=.4*Ly
    # Verify exact agreement with the independently saved daily movie diagnostic.
    movie=src/f'frame_{round(t/meta["output_seconds"]):05d}.jld2'
    with h5py.File(movie) as h:
        zmovie=np.asarray(h['video_heights']);qm=np.asarray(h['vorticity_slices']).transpose(2,1,0)
    iz=[int(np.argmin(abs(z-zz))) for zz in zmovie]
    movie_error=float(np.max(abs(q[:,:,iz]-qm)));assert movie_error<2e-10
    # Moments/PDFs use equal area and equal snapshot weights at each physical height.
    per=[];perface=[];perhist=[]
    for region in [mask,wall]:
        r=a[:,region,:];m1=r.mean((0,1));m2=(r*r).mean((0,1));m3=(r*r*r).mean((0,1))
        per.append(np.array([m1,m2,m3,r.min((0,1)),r.max((0,1)),(r< -1).mean((0,1)),(r>1).mean((0,1)),
                            (r<-.5).mean((0,1)),(r>.5).mean((0,1))]))
    rf=native[:,maskf,:];perface=np.array([rf.min((0,1)),rf.max((0,1)),(rf< -1).mean((0,1)),(rf>1).mean((0,1))])
    for k in range(nz):perhist.append(np.histogram(a[:,mask,k],bins=bins)[0])
    mom.append(per);face.append(perface);hists.append(perhist)
    weights=rho/rho.sum();U=np.sum(u*weights,axis=2);V=np.sum(v*weights,axis=2)
    up=u-U[:,:,None];vp=v-V[:,:,None]
    ux,uy,vx,vy=dx(U,hx),dy(U,hy),dx(V,hx),dy(V,hy,-1)
    # Positive P transfers horizontal depth-mean kinetic energy into deviations.
    production=-(up*up*ux[:,:,None]+up*vp*(uy+vx)[:,:,None]+vp*vp*vy[:,:,None])
    thbar=theta.mean((0,1));buoyproxy=meta['gravity']*w*(theta-thbar)/thbar
    ke3=(up*up+vp*vp+w*w)/2
    den=nx*int(mask.sum())*rho.sum()
    weighted=lambda arr:float(np.sum(arr[:,mask,:]*rho)/den)
    aa=a[:,mask,:];pp=production[:,mask,:];bb=buoyproxy[:,mask,:];kk=ke3[:,mask,:]
    conditional={}
    for name,sel in [('zeta_lt_minus_f',aa< -1),('zeta_gt_f',aa>1),('anticyclonic',aa<0),('cyclonic',aa>0)]:
        wt=np.broadcast_to(rho,sel.shape)*sel;mass=float(wt.sum())
        conditional[name]={'mass_fraction':mass/den,'shear_production':float(np.sum(wt*pp)/mass) if mass else None,
                           'buoyancy_proxy':float(np.sum(wt*bb)/mass) if mass else None,
                           'deviation_ke':float(np.sum(wt*kk)/mass) if mass else None}
    energy.append({'day':t/86400,'deviation_ke_m2_s2':weighted(ke3),'barotropic_to_deviations_m2_s3':weighted(production),
                   'buoyancy_production_proxy_m2_s3':weighted(buoyproxy),'conditional':conditional})
    for k in range(nz):
        for ij in candidates(u[:,:,k],v[:,:,k],q[:,:,k],f,hx,hy,mask,limit=3):
            r=radial(u[:,:,k],v[:,:,k],q[:,:,k],ij,f,hx,hy)
            c={'day':t/86400,'z_m':float(z[k]),'i':ij[0],'j':ij[1],
               'x_m':(ij[0]+.5)*hx,'y_m':float(y[ij[1]]),'core_zeta_over_f':float(a[ij[0],ij[1],k]),
               'coherence':r['coherence'],'min_Phi_over_f2_resolved':r['min_Phi_over_f2_resolved'],
               'min_Phi_circulation_over_f2_resolved':r['min_Phi_circulation_over_f2_resolved'],
               'resolved_negative_both':r['resolved_negative_both'],'axisymmetric_screen_pass':r['axisymmetric_screen_pass']}
            allcores.append(c);profiles.append(r)
    # A large-scale vortex view at the last time, with both signs retained.
    if p==files[-1]:
        Z=np.sum(q*weights,axis=2);lsvs=[]
        for sign in [-1,1]:
            inds=np.argwhere(np.broadcast_to(mask[None,:],Z.shape))
            ij=tuple(inds[np.argmax(sign*Z[:,mask].ravel())])
            r=radial(U,V,Z,ij,f,hx,hy)
            lsvs.append({'sign':sign,'x_m':(ij[0]+.5)*hx,'y_m':float(y[ij[1]]),'core_Z_over_f':float(Z[ij]/f),'profile':r})
        np.savez_compressed(out/'lsv_maps.npz',Z_over_f=Z/f,U=U,V=V,x=(np.arange(nx)+.5)*hx,y=y)
        save(out/'lsv_profiles.json',lsvs)
    meanN2=meta['gravity']*np.gradient(thbar,z)/thbar
    extremes.append({'day':t/86400,'min_zeta_over_f_interior':float(a[:,mask,:].min()),
        'max_zeta_over_f_interior':float(a[:,mask,:].max()),'minimum_native_corner_over_f':float(native[:,maskf,:].min()),
        'movie_curl_max_difference_s_inv':movie_error,'mean_N2_profile':meanN2})
    manifest.append({'file':p.name,'bytes':p.stat().st_size,'t_seconds':t,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    print('FULL_FIELD',p.name,'cores',len(allcores),'min',extremes[-1]['min_zeta_over_f_interior'],flush=True)

mom=np.asarray(mom);face=np.asarray(face);hists=np.asarray(hists)
assert np.allclose(np.diff([x['t_seconds'] for x in manifest]),86400)
aggregate=mom.mean(0);m1,m2,m3=aggregate[0,:3];raw=m3/np.maximum(m2,1e-30)**1.5
central=(m3-3*m1*m2+2*m1**3)/np.maximum(m2-m1*m1,1e-30)**1.5
pdf=hists.sum(0)/(len(files)*nx*int(mask.sum())*np.diff(bins)[None,:])
np.savez_compressed(out/'height_statistics.npz',z=z,bins=bins,pdf=pdf,moments=mom,native_face=face,raw_skewness=raw,central_skewness=central)
selected=sorted(range(len(allcores)),key=lambda i:(not allcores[i]['axisymmetric_screen_pass'],allcores[i]['core_zeta_over_f']))[:6]
save(out/'selected_vortex_profiles.json',[dict(allcores[i],profile=profiles[i]) for i in selected])
with (out/'vortex_candidates.csv').open('w') as fp:
    if allcores:
        wr=csv.DictWriter(fp,fieldnames=list(allcores[0]));wr.writeheader();wr.writerows(allcores)

# Three-hourly tracking only where full horizontal velocities are actually saved (~20 km).
tracks={};previous=[];nextid=0;timeline=[];rejected_links=0
for p in frames:
    with h5py.File(p) as h:
        u,v,q=[np.asarray(h[k]).T for k in ['u','v','zeta']];t=float(h['t'][()]);U=np.asarray(h['ubar']).T;V=np.asarray(h['vbar']).T
    cs=[]
    for sign in [-1,1]:
        for ij in candidates(u,v,q,f,hx,hy,mask,sign=sign,limit=8):
            r=radial(u,v,q,ij,f,hx,hy)
            cs.append({'time':t,'day':t/86400,'sign':sign,'i':ij[0],'j':ij[1],
                'x':(ij[0]+.5)*hx,'y':(ij[1]+.5)*hy,'u':float(U[ij]),'v':float(V[ij]),
                'q_over_f':float(q[ij]/f),'coherence':r['coherence'],
                'Phi_min':r['min_Phi_over_f2_resolved'],'screen':r['axisymmetric_screen_pass']})
    assigned=set()
    if previous and cs:
        cost=np.full((len(previous),len(cs)),1e6)
        for ia,a in enumerate(previous):
            dt=t-a['time'];xp=(a['x']+a['u']*dt)%Lx;yp=a['y']+a['v']*dt
            for ib,b in enumerate(cs):
                ddx=(b['x']-xp+Lx/2)%Lx-Lx/2;dist=np.hypot(ddx,b['y']-yp)
                ratio=abs(b['q_over_f']/a['q_over_f'])
                if a['sign']==b['sign'] and dist<=4*max(hx,hy) and .25<=ratio<=4:cost[ia,ib]=dist
        aa,bb=linear_sum_assignment(cost)
        for ia,ib in zip(aa,bb):
            if cost[ia,ib]<1e6:cs[ib]['id']=previous[ia]['id'];assigned.add(ib)
            else:rejected_links+=1
    for j,c in enumerate(cs):
        if j not in assigned:c['id']=nextid;nextid+=1
        tracks.setdefault(c['id'],[]).append(c)
    previous=cs
    timeline.append({'day':t/86400,'q_min':float(q[:,mask].min()/f),'q_max':float(q[:,mask].max()/f),
                     'fraction_lt_minus_f':float(np.mean(q[:,mask]/f< -1)),
                     'candidate_anticyclones':sum(c['sign']<0 for c in cs),'candidate_cyclones':sum(c['sign']>0 for c in cs),
                     'screened_anticyclones':sum(c['sign']<0 and c['screen'] for c in cs)})
save(out/'tracking.json',{'tracks':tracks,'timeline':timeline,'tracking_height_m':float(zmovie[0]),
     'rule':'At most eight extrema per sign/frame; rotation-dominated horizontal gradient; |zeta/f|>=0.25. One-to-one matching within four cells of a depth-mean-velocity-advected center, strength ratio 0.25–4.',
     'limitation':'Tracks are tentative pattern associations, not material vortices. Endpoints can mean weakening, merging, rank truncation or failed matching. Three-hour sampling cannot resolve fast centrifugal growth.'})
save(out/'energy.json',energy)
bands={'bottom_0_10km':z<10000,'lower_10_60km':(z>=10000)&(z<60000),
       'upper_convecting_60_80km':(z>=60000)&(z<80000),'cap_80_120km':z>=80000}
bandreport={}
for name,sel in bands.items():
    meanm=np.mean(aggregate[0,:3][:,sel],axis=1)
    bandreport[name]={'heights_m':z[sel].tolist(),'mean_zeta_over_f':float(meanm[0]),
                     'raw_skewness':float(meanm[2]/max(meanm[1],1e-30)**1.5),
                     'fraction_below_minus_f':float(mom[:,0,5,sel].mean()),
                     'fraction_above_f':float(mom[:,0,6,sel].mean()),
                     'native_corner_fraction_below_minus_f':float(face[:,2,sel].mean()),
                     'native_corner_fraction_above_f':float(face[:,3,sel].mean()),
                     'minimum_zeta_over_f':float(mom[:,0,3,sel].min()),'maximum_zeta_over_f':float(mom[:,0,4,sel].max()),
                     'pdf_outside_minus4_plus4_fraction':float(1-np.sum(pdf[sel]*np.diff(bins),axis=1).mean())}
passing=[c for c in allcores if c['axisymmetric_screen_pass']]
report={'status':'COMPLETE','compute_location':'Unity','source_run':str(src),'days':[90,100],
    'grid':[nx,ny,nz],'dx_m':hx,'dy_m':hy,'dz_m':hz,'f_s_inv':f,'cap_height_m':meta['cap_height'],
    'full_3d_snapshot_count':len(files),'full_3d_cadence_hours':24,'plane_snapshot_count':len(frames),'plane_cadence_hours':3,
    'inertial_period_hours':2*np.pi/f/3600,'inverse_f_hours':1/f/3600,
    'horizontal_statistics_region':'Central 80% of channel width; all periodic x. Wall strips measured separately.',
    'bands':bandreport,'extremes_by_day':extremes,'vortex_candidates':len(allcores),
    'axisymmetric_screen_pass_count':len(passing),'screened_candidates':passing,
    'profile_screen':'At least two adjacent radial samples at r=4–10 cells with BOTH discriminants < -0.01 f², and axisymmetric azimuthal energy / translated total horizontal energy >=0.5. Three strongest rotation-dominated anticyclones per height/time; not exhaustive.',
    'tracked_patterns':len(tracks),'tracking_height_m':float(zmovie[0]),
    'plane_screened_anticyclone_observations':sum(x['screened_anticyclones'] for x in timeline),
    'source_manifest':manifest,'validation':validation,'elapsed_seconds':time.time()-started,
    'limitations':['This tests a mature days-90–100 window, not the birth of the LSV.',
     'No smoothing is used in the PDFs. Native-corner and centered-curl differences measure operator sensitivity, not resolution convergence.',
     'Axisymmetric inviscid barotropic discriminants are conditional screens in this anelastic, convecting, capped flow; viscosity and vertical structure are not an eigenmode calculation.',
     'Daily 3D snapshots cannot resolve instability growth on the inverse-f time. Three-hourly horizontal fields exist only at the lowest movie level; other heights lack matching velocity histories.',
     'Centered Reynolds-stress production and the thermal covariance are diagnostic proxies, not a closed perturbation-energy budget. Positive production alone does not identify centrifugal instability.',
     'A suppressed negative tail alone cannot distinguish instability from asymmetric generation or numerical regularization.',
     'No new simulation, parameter change, or added instability forcing was used.']}
save(out/'analysis.json',report)
print('CENTRIFUGAL_ANALYSIS_COMPLETE',json.dumps({'bands':bandreport,'candidates':len(allcores),'screen_pass':len(passing),'elapsed_seconds':report['elapsed_seconds']}),flush=True)

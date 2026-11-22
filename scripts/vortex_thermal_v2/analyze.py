"""Saved-field, read-only analysis. All field processing must run on Unity."""
import sys, json, csv, hashlib
from pathlib import Path
import numpy as np
import h5py
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
L=4.; DC=.78; CP=2.5

def periodic_components(mask):
    lab,n=label(mask); parent=np.arange(n+1)
    def root(i):
        while parent[i]!=i: parent[i]=parent[parent[i]];i=parent[i]
        return i
    for a,b in list(zip(lab[0],lab[-1]))+list(zip(lab[:,0],lab[:,-1])):
        if a and b: parent[root(a)]=root(b)
    lut=np.array([root(i) for i in range(n+1)])
    return lut[lab]

def core(zeta,sign,previous=None,threshold=.4):
    mag=sign*zeta; peak=mag.max()
    if peak<1e-9:return None,dict(valid=False)
    labs=periodic_components(mag>threshold*peak); ids=np.unique(labs);ids=ids[ids>0]
    candidates=[]
    for i in ids:
        m=labs==i
        overlap=0. if previous is None else np.sum(m&previous)/max(np.sum(m|previous),1)
        candidates.append((overlap,float(mag[m].sum()),m))
    choice=max(candidates,key=lambda q:(q[0],q[1])) if previous is not None else max(candidates,key=lambda q:q[1])
    ov,_,m=choice;n=zeta.shape[0];dx=L/n
    xc=(np.arange(n)+.5)*dx;yc=xc
    def center(weights):return float(np.mod(np.angle(np.sum(weights*np.exp(2j*np.pi*xc/L))),2*np.pi)*L/(2*np.pi))
    return m,dict(valid=True,x=center(m.sum(0)),y=center(m.sum(1)),area=float(m.sum()*dx*dx),circulation=float(zeta[m].sum()*dx*dx),peak=float(sign*mag[m].max()),previous_iou=float(ov),identity_ambiguous=bool(previous is not None and ov<.2))

def low_vorticity(u,v):
    n=u.shape[-1];dx=L/n
    z=(np.roll(v,-1,-1)-np.roll(v,1,-1)-np.roll(u,-1,-2)+np.roll(u,1,-2))/(2*dx)
    k=np.fft.fftfreq(n)*n;kx,ky=np.meshgrid(k,k);kr=np.hypot(kx,ky)
    return z,np.fft.ifft2(np.fft.fft2(z)*((kr>0)&(kr<=2))).real,kr

def checks():
    n=64;x=(np.arange(n)+.5)*L/n;xx,yy=np.meshgrid(x,x)
    z=np.cos(2*np.pi*xx/L)*np.cos(2*np.pi*yy/L)
    m,a=core(z,1);assert a['valid'] and a['circulation']>0
    mask=np.zeros((n,n),bool);mask[10:20,:3]=True;mask[10:20,-3:]=True
    assert len(np.unique(periodic_components(mask)[mask]))==1
    dx=L/n;u=-np.sin(2*np.pi*yy/L);v=np.sin(2*np.pi*xx/L)
    zz,zl,_=low_vorticity(u,v);want=(np.cos(2*np.pi*xx/L)+np.cos(2*np.pi*yy/L))*np.sin(2*np.pi/n)/dx
    assert np.max(abs(zz-want))<1e-12 and np.max(abs(zz-zl))<1e-12
    assert abs(zz.sum()*dx*dx)<1e-12
    print('PERIODIC_MASK_CURL_FILTER_CHECKS_PASS',flush=True)

def writecsv(path,rows):
    with path.open('w') as h:
        w=csv.DictWriter(h,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));w.writeheader();w.writerows(rows)

def main():
    checks();out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    roots={'omega3':['data/chan-eddies-v3-screen-2089476/case4','data/chan-eddies-v4-case4-2090026/parent','data/chan-rotation-response-v1-control-2093290/parent','data/omega3-long-v1-run-2096858/parent'], 'omega45':['data/chan-fixed45-v2-run-2093564/parent']}
    summary={}; provenance={}
    for name,segments in roots.items():
        omega=3 if name=='omega3' else 4.5;frames={};fields={}
        for seg in map(Path,segments):
            for p in seg.glob('frame_*.jld2'):
                if 'control-2093290' in str(p) and p.name=='frame_0026307.jld2':continue
                with h5py.File(p) as h:t=float(h['t'][()])
                if abs(t-round(t))<1e-8:frames[t]=p
            for p in (seg/'fields3d').glob('*.jld2'):
                with h5py.File(p) as h:t=float(h['t'][()])
                if abs(t-263.06682088232816)<1e-6:continue
                fields[t]=p
        previous={1:None,-1:None};tracks=[];profiles=[];convection=[];parity=[]
        provenance[name]={'segments':segments,'frames':len(frames),'fields3d':len(fields),'sample_times':sorted(fields)}
        for t,p in sorted(frames.items()):
            with h5py.File(p) as h:
                ub=h['ubar'][()].astype(float);vb=h['vbar'][()].astype(float);f=float(h['f'][()])
            zb,zl,kr=low_vorticity(ub,vb);regions={}
            for sign,region in [(1,'cyclone'),(-1,'anticyclone')]:
                m,rec=core(zl,sign,previous[sign]);previous[sign]=m;regions[region]=m
                tracks.append(dict(time=t,region=region,**rec))
                for th in [.25,.55]:
                    _,r=core(zl,sign,threshold=th)
                    tracks[-1].update({f'area_threshold_{th}':r.get('area',0),f'circulation_threshold_{th}':r.get('circulation',0)})
            if any(m is None for m in regions.values()):continue
            regions['background']=~(regions['cyclone']|regions['anticyclone'])
            if t not in fields:continue
            with h5py.File(fields[t]) as h:
                z=h['z'][()].astype(float);rho=h['density'][()].astype(float);T=h['temperature'][()].astype(float)
                uf=h['u_face'][()].astype(float);vf=h['v_face'][()].astype(float);wf=h['w_face'][()].astype(float)
            nz,ny,nx=rho.shape;assert nx==ny==256 and nz==80 and wf.shape==(nz+1,ny,nx)
            assert rho.min()>0 and T.min()>0 and all(np.isfinite(a).all() for a in [rho,T,uf,vf,wf])
            u=.5*(uf+np.roll(uf,-1,2));v=.5*(vf+np.roll(vf,-1,1));w=.5*(wf[:-1]+wf[1:])
            den=rho.sum(0);uavg=(rho*u).sum(0)/den;vavg=(rho*v).sum(0)/den
            err=max(float(abs(ub-uavg).max()),float(abs(vb-vavg).max()));assert err<2e-6,err;parity.append(err)
            Tmean=T.mean((1,2))[:,None,None];rmean=rho.mean((1,2))[:,None,None]
            wt=(rho*w).mean((1,2))[:,None,None]/rmean;Tt=(rho*T).mean((1,2))[:,None,None]/rmean
            Ta=T-Tmean;wp=w-wt;Tp=T-Tt
            heat=CP*rho*wp*Tp;heat_anomaly=CP*rho*w*Ta;raw_heat=CP*rho*w*T
            fluct=[]
            for cutoff in [2,3]:
                mask=kr>cutoff;ener=np.zeros_like(kr);variance=np.zeros(nz);hp=[]
                for vel in [u,v,w]:
                    vh=np.fft.fft2(vel,axes=(-2,-1));vh*=mask
                    vv=np.fft.ifft2(vh,axes=(-2,-1)).real
                    variance+=np.mean(vv**2,(1,2));ener+=np.mean(abs(vh[z<DC])**2,axis=0)
                    if cutoff==2:hp.append(vv)
                kmean=float(np.sum(ener*kr)/max(ener.sum(),1e-99));ell=L/kmean if kmean else 0
                rms=float(np.sqrt(variance[z<DC].mean()));ro=rms/(2*omega*ell) if ell else 0
                convection.append(dict(time=t,cutoff=cutoff,rms=rms,horizontal_wavelength=ell,rossby_2omega=ro,rossby_f=rms/(f*ell) if ell else 0))
                if cutoff==2:fluct=hp
            for region,m in regions.items():
                for k,zk in enumerate(z):
                    profiles.append(dict(time=t,region=region,z=float(zk),area_fraction=float(m.mean()),temperature_anomaly=float(Ta[k,m].mean()),w_rms=float(np.sqrt(np.mean(w[k,m]**2))),w_fluct_rms=float(np.sqrt(np.mean(fluct[2][k,m]**2))),convective_enthalpy_flux=float(heat[k,m].mean()),enthalpy_temperature_anomaly=float(heat_anomaly[k,m].mean()),raw_enthalpy_flux=float(raw_heat[k,m].mean()),mass_flux=float((rho[k]*w[k])[m].mean())))
            print(name,'FIELD_ANALYZED',t,flush=True)
        assert profiles and convection
        writecsv(out/f'{name}_vortex_tracks.csv',tracks);writecsv(out/f'{name}_conditional_profiles.csv',profiles);writecsv(out/f'{name}_convective_scales.csv',convection)
        times=sorted(fields);late=[r for r in profiles if r['time']>=max(times)-100];matched=[r for r in profiles if 400<=r['time']<=480]
        def means(rows):
            return {reg:{key:float(np.mean([r[key] for r in rows if r['region']==reg and r['z']<DC])) for key in ['temperature_anomaly','w_fluct_rms','convective_enthalpy_flux','mass_flux']} for reg in regions}
        summary[name]={'last_time':max(frames),'last_3d_time':max(times),'field_velocity_parity_max':max(parity),'last_100_bulk_means':means(late),'matched_400_480_bulk_means':means(matched),'last_100_convective_scales':{key:float(np.mean([r[key] for r in convection if r['time']>=max(times)-100 and r['cutoff']==2])) for key in ['rms','horizontal_wavelength','rossby_2omega','rossby_f']},'ambiguous_identity_samples':sum(r.get('identity_ambiguous',False) for r in tracks)}
        fig,axs=plt.subplots(2,3,figsize=(12,7))
        for reg,col in [('cyclone','firebrick'),('anticyclone','royalblue'),('background','gray')]:
            rr=[r for r in tracks if r['region']==reg and r.get('valid',False)]
            if rr:
                axs[0,0].plot([r['time'] for r in rr],[r['area'] for r in rr],c=col,label=reg)
                axs[0,1].plot([r['time'] for r in rr],[r['circulation'] for r in rr],c=col)
            for j,key in enumerate(['temperature_anomaly','w_fluct_rms','convective_enthalpy_flux']):
                prof=[np.mean([r[key] for r in late if r['region']==reg and r['z']==zk]) for zk in z]
                axs[1,j].plot(prof,z,c=col,label=reg)
                axs[1,j].axhline(DC,color='k',ls=':',lw=.6);axs[1,j].set_xlabel(key);axs[1,j].set_ylabel('z/H')
        for cut in [2,3]:
            rr=[r for r in convection if r['cutoff']==cut];axs[0,2].plot([r['time'] for r in rr],[r['rossby_2omega'] for r in rr],label=f'|k|>{cut}')
        axs[0,0].set_ylabel('Core area');axs[0,1].set_ylabel('Signed filtered circulation');axs[0,2].set_ylabel('Fluctuation Rossby number');axs[0,0].legend();axs[0,2].legend()
        for ax in axs[0]:ax.set_xlabel('Dimensionless time')
        fig.suptitle(f'{name}: tracked filtered cores; lower row averages last 100 saved time units');fig.tight_layout();fig.savefig(out/f'{name}_thermal.png',dpi=140);plt.close(fig)
    methods={'core':'Curl of density-weighted depth-mean velocity filtered to 0<|k|<=2; periodic connected regions above 40% of each signed peak; overlap tracking, ambiguous below IoU=.2. Threshold sensitivity .25,.55 retained. Core columns used at every height, not independently tracked slice vortices.','flux':'Cp*rho*(w-Favre_plane_mean(w))*(T-Favre_plane_mean(T)), Cp=2.5. Also save raw enthalpy, plane-T-anomaly flux and mass flux to expose mean advection.','rossby':'Remove horizontal |k|<=2 (also <=3 sensitivity) separately at each depth for u,v,w. RMS in z<.78. Horizontal wavelength=4/energy-weighted integer k. Ro=uprime/(2Omega*wavelength), also f convention. Definition explicitly operational, not assumed equal to any paper threshold.','limitations':['One realization per rotation; finite-amplitude initial eddies.','Energy budget is unclosed; conditional correlations cannot establish causality or a closed maintenance budget.','Core membership threshold and identity ambiguity must be considered.','3D samples every20; tracks every1; exclude known invalid restart sample.']}
    (out/'summary.json').write_text(json.dumps({'cases':summary,'methods':methods,'provenance':provenance},indent=2)+'\n')
    (out/'ANALYSIS_PASS.json').write_text(json.dumps({'status':'PASS','analysis_on':'Unity','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':list(summary)},indent=2)+'\n')
    print('ANALYSIS_COMPLETE',json.dumps(summary),flush=True)
if __name__=='__main__':main()

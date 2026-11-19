"""Read-only analysis of completed runs. Execute on a Unity compute node only.

Centered reconstructions are explicitly NOT a native discrete tendency budget.
"""
import csv, hashlib, json, os, sys, time, tomllib
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = {'none':1990851, 'cyclonic':1990854, 'anticyclonic':1990856}
def dx(a,h): return (np.roll(a,-1,axis=-1)-np.roll(a,1,axis=-1))/(2*h)
def dy(a,h,odd=False):
    sign=-1 if odd else 1
    p=np.concatenate((sign*a[...,0:1,:],a,sign*a[...,-1:,:]),axis=-2)
    return (p[...,2:,:]-p[...,:-2,:])/(2*h)
def nonzonal(a): return a-a.mean(axis=-1,keepdims=True)
def skew(a): return float(np.mean(a**3)/max(np.mean(a*a)**1.5,1e-50))
def rms(a): return float(np.sqrt(np.mean(a*a)))
def proj(a,Lx,Ly,cutoff,odd=False):
    ny,nx=a.shape
    ae=np.concatenate((a,(-1 if odd else 1)*a[::-1]),axis=0)
    ah=np.fft.fft2(ae,norm='forward')
    kx=2*np.pi*np.fft.fftfreq(nx,Lx/nx)[None,:]
    ky=2*np.pi*np.fft.fftfreq(2*ny,Ly/ny)[:,None]
    mask=(np.hypot(kx,ky)<2*np.pi/cutoff)&(kx!=0)
    assert np.isclose(np.sum(abs(ah)**2),np.mean(ae*ae),rtol=1e-10,atol=1e-15)
    return np.fft.ifft2(ah*mask,norm='forward').real[:ny]
def writecsv(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader();w.writerows(rows)
def test_operators():
    nx,ny=64,40; x=(np.arange(nx)+.5)*2*np.pi/nx
    a=np.broadcast_to(np.sin(3*x),(ny,nx)); expected=3*np.cos(3*x)*np.sin(3*2*np.pi/nx)/(3*2*np.pi/nx)
    assert np.max(abs(dx(a,2*np.pi/nx)-expected))<1e-12
    assert np.max(abs(dy(np.ones_like(a),1)))==0
    u=np.broadcast_to(x,(ny,nx));v=np.broadcast_to(np.arange(ny)[:,None],(ny,nx))
    assert np.max(abs(dx(v,1)-dy(u,1)))==0
    q=proj(a,2*np.pi,2*np.pi,1)
    assert np.max(abs(q-a))<1e-12
    return 'PASS analytic derivative, zero curl, constant-wall extension, spectral projection and Parseval'

def main():
    assert 'SLURM_JOB_ID' in os.environ, 'Run analysis in a Unity allocation'
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=False)
    started=time.time();validation={'operators':test_operators()}
    profiles=[]; layers=[]; history=[]; summary={}; inventory=[]
    for case,job in CASES.items():
        root=Path(f'data/strong256-direct-v1-{job}'); seg=root/'day1_50'
        m=tomllib.loads((seg/'metadata.toml').read_text())
        assert m['status']=='COMPLETE' and m['simulated_seconds']==50*86400 and m['beta']==0
        nx,ny,nz=m['nx'],m['ny'],m['nz']; Lx,Ly,H=m['lx'],m['ly'],m['height'];hx,hy,hz=Lx/nx,Ly/ny,H/nz
        f=m['f']; cut=2*m['turbulence_wavelength']; core=slice(int(np.ceil(.1*ny))-1,int(np.floor(.9*ny)))
        initial=root/'day0_1/frame_00000.jld2'
        with h5py.File(initial) as d: initial_jet=d['ubar'][:].mean(axis=-1)
        shear_sign=m['shear_rossby']
        nframes=0
        for p in sorted(seg.glob('frame_*.jld2')):
            with h5py.File(p) as d:
                day=float(d['t'][()])/86400
                if day<40: continue
                u=d['ubar'][:];v=d['vbar'][:]
                Rxx=d['Rxx'][:];Rxy=d['Rxy'][:];Ryy=d['Ryy'][:]
                ul=proj(nonzonal(u),Lx,Ly,cut);vl=proj(nonzonal(v),Lx,Ly,cut,True)
                fx=-dx(Rxx,hx)-dy(Rxy,hy,True);fy=-dx(Rxy,hx)-dy(Ryy,hy)
                ax=-u*dx(u,hx)-v*dy(u,hy);ay=-u*dx(v,hx)-v*dy(v,hy,True)
                jet=u.mean(axis=-1); j=jet[core]; j0=initial_jet[core]
                history.append(dict(case=case,day=day,large_nonzonal_ke_m2_s2=float(.5*np.mean(ul*ul+vl*vl)),
                    resolved_horizontal_baroclinic_stress_work_m2_s3=float(np.mean(ul*fx+vl*fy)),
                    resolved_barotropic_advection_work_m2_s3=float(np.mean(ul*ax+vl*ay)),
                    zonal_rms_m_s=rms(j),initial_zonal_rms_m_s=rms(j0),
                    zonal_projection_on_initial=float(np.dot(j,j0)/np.dot(j0,j0)) if shear_sign else None,
                    zonal_vorticity_core_mean_over_f=float(np.mean(-np.gradient(jet,hy)[core])/f)))
                nframes+=1
        assert nframes==81,(case,nframes)
        max_movie_error=0.; max_column_massdiv_ratio=0.;dates=[];grid_metrics=[]
        for p in sorted((seg/'fields3d').glob('frame_*.jld2')):
            with h5py.File(p) as d:
                day=float(d['t'][()])/86400
                if day<40: continue
                dates.append(day);rho=d['density'][:];z=d['z'][:]
                uf=d['u_face'][:].astype(float);vf=d['v_face'][:].astype(float);wf=d['w_face'][:].astype(float)
                th=d['theta'][:].astype(float)
            assert uf.shape==(nz,ny,nx) and vf.shape==(nz,ny+1,nx) and wf.shape==(nz+1,ny,nx)
            assert np.max(abs(vf[:,[0,-1]]))==0 and np.max(abs(wf[[0,-1]]))==0
            u=.5*(uf+np.roll(uf,-1,axis=2)); v=.5*(vf[:,1:]+vf[:,:-1]);w=.5*(wf[1:]+wf[:-1])
            div=(np.roll(uf,-1,axis=2)-uf)/hx+(vf[:,1:]-vf[:,:-1])/hy
            divc=dx(u,hx)+dy(v,hy,True)
            zz=dx(v,hx)-dy(u,hy);zp=nonzonal(zz)
            stretch=-(f+zz)*divc
            tilt=np.gradient(u,hz,axis=0)*dy(w,hy)-np.gradient(v,hz,axis=0)*dx(w,hx)
            wc=w[:,core];dc=divc[:,core];zc=zp[:,core]
            wr=np.sqrt(np.mean(wc*wc,axis=(1,2)))
            col=np.sum(rho[:,None,None]*div,axis=0)
            denom=np.sum(rho*np.sqrt(np.mean(div**2,axis=(1,2))))
            max_column_massdiv_ratio=max(max_column_massdiv_ratio,rms(col)/max(denom,1e-50))
            # Reconstruct the actual movie operator independently from saved face velocities.
            with h5py.File(seg/p.name) as d:
                idx=[int(np.argmin(abs(z-h))) for h in d['video_heights'][:]]
                saved=d['vorticity_slices'][:].astype(float)
                max_movie_error=max(max_movie_error,float(np.max(abs(zz[idx]-saved))/max(np.max(abs(saved)),1e-30)))
                assert max_movie_error<2e-5, 'Saved-field curl mismatch'
            for k,zk in enumerate(z):
                up=wc[k]>wr[k];down=wc[k]<-wr[k]
                cond=lambda a,mask:float(np.mean(a[mask])) if np.any(mask) else None
                theta_anom=th[k,core]-th[k,core].mean()
                profiles.append(dict(case=case,day=day,z_km=float(zk/1000),
                    nonzonal_zeta_rms_over_f=rms(zc[k])/abs(f),signed_nonzonal_skew=skew(zc[k])*np.sign(f),
                    w_rms=wr[k],up_area=float(up.mean()),down_area=float(down.mean()),
                    up_divergence_over_f=cond(dc[k]/abs(f),up),down_divergence_over_f=cond(dc[k]/abs(f),down),
                    up_zeta_prime_over_f=cond(zc[k]/f,up),down_zeta_prime_over_f=cond(zc[k]/f,down),
                    up_stretch_over_f2=cond(stretch[k,core]/f**2,up),up_tilt_over_f2=cond(tilt[k,core]/f**2,up),
                    theta_anomaly_up_K=cond(theta_anom,up),
                    stretching_M3_rate=float(3*np.mean(zc[k]**2*nonzonal(stretch[k])[core]))*np.sign(f),
                    tilting_M3_rate=float(3*np.mean(zc[k]**2*nonzonal(tilt[k])[core]))*np.sign(f)))
            for label,mask in [('lower',(z>=10000)&(z<=70000)),('interface',(z>=75000)&(z<=90000)),('upper',(z>=90000)&(z<=105000)),('full',np.ones(nz,dtype=bool))]:
                weights=rho[mask]/rho[mask].sum();avg=lambda a:np.sum(weights[:,None,None]*a[mask],axis=0)
                za=avg(zz); pp=nonzonal(za); up=avg(u);vp=avg(v)
                ul=proj(nonzonal(up),Lx,Ly,cut);vl=proj(nonzonal(vp),Lx,Ly,cut,True)
                layers.append(dict(case=case,day=day,layer=label,actual_z_min_km=float(z[mask].min()/1000),actual_z_max_km=float(z[mask].max()/1000),
                    signed_skew_nonzonal=skew(pp[:,core])*np.sign(f),rms_zeta_nonzonal_over_f=rms(pp[:,core])/abs(f),
                    rms_zeta_total_over_f=rms(za[:,core])/abs(f),nonzonal_ke=float(.5*np.mean(nonzonal(up)**2+nonzonal(vp)**2)),
                    large_nonzonal_ke=float(.5*np.mean(ul**2+vl**2))))
            for height in (20000,60000,80000,90000,100000):
                k=int(np.argmin(abs(z-height)));a=nonzonal(w[k]);ae=np.concatenate((a,a[::-1]),axis=0);A=np.fft.fft2(ae,norm='forward')
                kx=2*np.pi*np.fft.fftfreq(nx,hx)[None,:];ky=2*np.pi*np.fft.fftfreq(2*ny,hy)[:,None];K=np.hypot(kx,ky);power=abs(A)**2
                grid_metrics.append(dict(day=day,z_km=float(z[k]/1000),w_variance_wavelength_below_4dx=float(power[K>2*np.pi/(4*hx)].sum()/max(power.sum(),1e-50))))
            inventory.append(dict(case=case,path=str(p),bytes=p.stat().st_size,time_days=day))
            print(f'READ {case} day={day:g}',flush=True)
        assert dates==list(range(40,51)),dates
        hh=[r for r in history if r['case']==case];ll=[r for r in layers if r['case']==case]
        summary[case]=dict(initial_shear_rossby=shear_sign,initial_jet_rms_m_s=hh[0]['initial_zonal_rms_m_s'],
            late_jet_rms_m_s=float(np.mean([r['zonal_rms_m_s'] for r in hh])),
            late_jet_projection_on_initial=float(np.mean([r['zonal_projection_on_initial'] for r in hh])) if shear_sign else None,
            late_horizontal_stress_work=float(np.mean([r['resolved_horizontal_baroclinic_stress_work_m2_s3'] for r in hh])),
            late_barotropic_advection_work=float(np.mean([r['resolved_barotropic_advection_work_m2_s3'] for r in hh])),
            layers={label:{key:float(np.mean([r[key] for r in ll if r['layer']==label])) for key in ('signed_skew_nonzonal','rms_zeta_nonzonal_over_f','large_nonzonal_ke')} for label in ('lower','interface','upper','full')},
            vertical_velocity_resolution=grid_metrics,
            checks=dict(movie_curl_max_relative_difference=max_movie_error,column_mass_divergence_relative=max_column_massdiv_ratio))
    writecsv(out/'profiles.csv',profiles);writecsv(out/'layers.csv',layers);writecsv(out/'history.csv',history)
    report=dict(status='COMPLETE',window_days=[40,50],daily_volume_states_per_case=11,layer_history_states_per_case=81,
        compute_location='Unity',slurm_job=os.environ['SLURM_JOB_ID'],analysis_seconds=time.time()-started,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),validation=validation,cases=summary,inputs=inventory,
        limitations=['Single initial seed; time samples are correlated.',
        'Up/down events use w > RMS(w) and w < -RMS(w), with RMS on fixed interior y cells. Conditional averages do not track plume trajectories.',
        'Centered stretching and tilting are instantaneous reconstructed terms, not native WENO/RK tendencies. No closed layerwise causal budget is claimed.',
        'Daily 3D fields cannot resolve individual plume growth or rapid instability. Three-hourly barotropic fields do not supply missing native terms.',
        'Large-scale projection uses wavelengths >2400 km, nonzonal kx != 0, and free-slip parity across channel walls. Stress work is a resolved reconstruction, not proof of an inverse cascade.',
        'Full-column horizontal baroclinic stress work excludes decomposition corrections from numerical operators; no complete kinetic-energy budget is claimed.',
        'Mass-continuity check uses native horizontal divergence and saved Float32 face velocities. It is not a vorticity budget.',
        'Near-grid vertical-velocity power is a warning about resolved plume physics, not a measured physical plume size.'])
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    fig,axs=plt.subplots(1,3,figsize=(13,5),layout='constrained',sharey=True)
    for case in CASES:
        rr=[r for r in profiles if r['case']==case];zz=sorted(set(r['z_km'] for r in rr))
        for ax,key in zip(axs,('signed_nonzonal_skew','nonzonal_zeta_rms_over_f','up_divergence_over_f')):
            vals=[np.mean([r[key] for r in rr if r['z_km']==z]) for z in zz]
            ax.plot(vals,zz,label=case);ax.axhline(80,color='.6',ls='--');ax.axvline(0,color='.7',lw=.5)
    for ax,label in zip(axs,('Nonzonal vorticity skewness','RMS nonzonal vorticity / |f|','Divergence / |f| in upflows')):ax.set_xlabel(label)
    axs[0].set_ylabel('Height above bottom [km]');axs[2].legend();fig.suptitle('Existing 256² Breeze shear cases · daily fields, Earth days 40–50\nInterior statistics; dashed line is the original stable-layer base')
    fig.savefig(out/'vertical_structure.png',dpi=150);plt.close(fig)
    print('ANALYSIS_COMPLETE '+str(out),flush=True)
if __name__=='__main__':main()

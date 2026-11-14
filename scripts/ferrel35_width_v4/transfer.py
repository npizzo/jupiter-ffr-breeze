"""Resolved depth-mean transfers in a free-slip channel; run on Unity.

Not an exact discrete WENO kinetic-energy budget. Negative cumulative advective
flux means upscale transfer. Stress work is a separate exchange term, not a flux.
"""
import argparse,csv,json,tomllib,hashlib
from pathlib import Path
import numpy as np
import h5py


def diagnose(u,v,rxx,rxy,ryy,lx,ly,wavelength):
    ny,nx=u.shape
    assert all(a.shape==(ny,nx) and np.isfinite(a).all() for a in (u,v,rxx,rxy,ryy))
    # Reflect cell-centered fields about both free-slip walls. The extension is
    # periodic over 2 Ly; means there equal the physical-channel spatial means.
    extend=lambda a,parity:np.concatenate((a,parity*a[::-1,:]),axis=0)
    fft=lambda a:np.fft.fft2(a,norm='forward')
    inv=lambda a:np.fft.ifft2(a,norm='forward').real
    U0,V0=fft(extend(u,1)),fft(extend(v,-1))
    kx=2*np.pi*np.fft.fftfreq(nx,lx/nx)[None,:]
    ky=2*np.pi*np.fft.fftfreq(2*ny,ly/ny)[:,None]
    K2=kx*kx+ky*ky;K=np.sqrt(K2);den=np.where(K2>0,K2,1)
    mask=(np.abs(np.fft.fftfreq(nx)*nx)[None,:]<nx/3)&(np.abs(np.fft.fftfreq(2*ny)*2*ny)[:,None]<2*ny/3)
    U,V=U0*mask,V0*mask
    filtered_energy=.5*np.sum(abs(U)**2+abs(V)**2)
    divergent=(kx*U+ky*V)/den
    U,V=U-kx*divergent,V-ky*divergent
    energy=.5*(abs(U)**2+abs(V)**2)
    ur,vr=inv(U),inv(V)
    Nu=-fft(ur*inv(1j*kx*U)+vr*inv(1j*ky*U))*mask
    Nv=-fft(ur*inv(1j*kx*V)+vr*inv(1j*ky*V))*mask
    # Projecting Nu,Nv is unnecessary for their inner product with solenoidal U,V.
    transfer=(np.conj(U)*Nu+np.conj(V)*Nv).real
    X,Y,Z=fft(extend(rxx,1))*mask,fft(extend(rxy,-1))*mask,fft(extend(ryy,1))*mask
    stress=(np.conj(U)*(-1j*kx*X-1j*ky*Y)+np.conj(V)*(-1j*kx*Y-1j*ky*Z)).real
    original_energy=.5*np.mean(extend(u,1)**2+extend(v,-1)**2)
    assert np.isclose(original_energy,.5*np.sum(abs(U0)**2+abs(V0)**2),rtol=1e-12,atol=1e-14)
    # Normalize near-zero analytic transfers by a physical advective rate;
    # dividing roundoff by roundoff gives a misleading O(1) relative error.
    rate_scale=np.mean(ur*ur+vr*vr)**1.5*np.sqrt(np.sum(K2*energy)/max(energy.sum(),1e-30))
    closure=abs(transfer.sum())/max(np.abs(transfer).sum(),rate_scale,1e-30)
    assert closure<1e-10,closure
    dk=min(2*np.pi/lx,np.pi/ly)
    bins=np.floor(K/dk+1e-10).astype(int)
    shells=lambda a:np.bincount(bins.ravel(),weights=a.ravel())
    E,T,S=shells(energy),shells(transfer),shells(stress)
    nonzero=K>0;k0=2*np.pi/wavelength
    large=nonzero&(K<.75*k0)
    metrics={
        'depth_mean_ke_raw_m2_s2':float(original_energy),
        'depth_mean_ke_analyzed_m2_s2':float(energy.sum()),
        'retained_energy_fraction':float(energy.sum()/max(original_energy,1e-30)),
        'projection_removed_fraction_of_filtered_energy':float(max(0,filtered_energy-energy.sum())/max(filtered_energy,1e-30)),
        'advective_transfer_closure_relative':float(closure),
        'large_scale_ke_m2_s2':float(energy[large].sum()),
        'large_scale_fraction_excluding_mean':float(energy[large].sum()/max(energy[nonzero].sum(),1e-30)),
        'flux_below_seed_band_m2_s3':float(-transfer[large].sum()),
        'stress_work_below_seed_band_m2_s3':float(stress[large].sum()),
        'stress_work_total_m2_s3':float(stress.sum()),
        'energy_weighted_inverse_wavenumber_m':float((energy[nonzero]/K[nonzero]).sum()/max(energy[nonzero].sum(),1e-30)),
        'zonal_ke_fraction':float(energy[np.broadcast_to(kx==0,energy.shape)&nonzero].sum()/max(energy[nonzero].sum(),1e-30))}
    return metrics,dict(k_upper=(np.arange(len(E))+1)*dk,energy_shell=E,advective_transfer_shell=T,advective_flux=-np.cumsum(T),stress_work_shell=S)


def self_test():
    ny=nx=64;lx=2*np.pi;ly=np.pi
    x=(np.arange(nx)+.5)*lx/nx;y=(np.arange(ny)+.5)*ly/ny
    X,Y=np.meshgrid(x,y)
    # A single Laplacian eigenmode has zero solenoidal nonlinear transfer.
    u=np.cos(Y)*np.cos(2*X);v=2*np.sin(Y)*np.sin(2*X);z=np.zeros_like(u)
    m,s=diagnose(u,v,z,z,z,lx,ly,2*np.pi/4)
    assert np.isclose(m['depth_mean_ke_analyzed_m2_s2'],.625,atol=1e-12)
    assert np.max(abs(s['advective_flux']))<1e-12
    # Unequal-scale interacting triads: closure, nonzero transfer, cubic scaling,
    # reversal of transfer under velocity reversal, and quadratic energy scaling.
    u=u+.7*3*np.cos(3*Y)*np.sin(X)+.4*2*np.cos(2*Y)*np.cos(3*X+.3)
    v=v-.7*np.sin(3*Y)*np.cos(X)+.4*3*np.sin(2*Y)*np.sin(3*X+.3)
    m,s=diagnose(u,v,z,z,z,lx,ly,2*np.pi/4)
    assert m['projection_removed_fraction_of_filtered_energy']<1e-12
    assert np.max(abs(s['advective_flux']))>1e-5
    _,neg=diagnose(-u,-v,z,z,z,lx,ly,2*np.pi/4)
    m2,twice=diagnose(2*u,2*v,z,z,z,lx,ly,2*np.pi/4)
    np.testing.assert_allclose(neg['advective_flux'],-s['advective_flux'],atol=1e-12)
    np.testing.assert_allclose(twice['advective_flux'],8*s['advective_flux'],atol=1e-11)
    np.testing.assert_allclose(twice['energy_shell'],4*s['energy_shell'],atol=1e-12)
    # Positive eddy-viscosity stresses remove energy from the resolved mean.
    ue= np.concatenate((u,u[::-1]),axis=0);ve=np.concatenate((v,-v[::-1]),axis=0)
    kx=2*np.pi*np.fft.fftfreq(nx,lx/nx)[None,:];ky=2*np.pi*np.fft.fftfreq(2*ny,ly/ny)[:,None]
    dx=lambda a:np.fft.ifft2(1j*kx*np.fft.fft2(a)).real
    dy=lambda a:np.fft.ifft2(1j*ky*np.fft.fft2(a)).real
    rxx=-2*dx(ue)[:ny];rxy=-(dy(ue)+dx(ve))[:ny];ryy=-2*dy(ve)[:ny]
    ms,_=diagnose(u,v,rxx,rxy,ryy,lx,ly,2*np.pi/4)
    assert ms['stress_work_total_m2_s3']<0
    expected=np.mean(dx(ue)*np.concatenate((rxx,rxx[::-1]))+(dy(ue)+dx(ve))*np.concatenate((rxy,-rxy[::-1]))+dy(ve)*np.concatenate((ryy,ryy[::-1])))
    assert np.isclose(ms['stress_work_total_m2_s3'],expected,rtol=1e-12)
    print('PASS: channel parity, known-mode energy, Parseval, conservative nonlinear transfer, transfer reversal/scaling, stress-work sign and integration by parts.',flush=True)


def analyze(folder):
    m=tomllib.loads((folder/'metadata.toml').read_text())
    assert m['status']=='COMPLETE'
    rows=[];spectra=[]
    for path in sorted(folder.glob('frame_*.jld2')):
        with h5py.File(path) as d:
            arrays=[d[key][:] for key in ['ubar','vbar','Rxx','Rxy','Ryy']]
            assert arrays[0].shape==(m['ny'],m['nx'])
            metrics,s=diagnose(*arrays,m['lx'],m['ly'],m['turbulence_wavelength'])
            metrics={'time_seconds':float(d['t'][()]),**metrics}
        rows.append(metrics);spectra.append(s)
    assert rows and all(b['time_seconds']>a['time_seconds'] for a,b in zip(rows,rows[1:]))
    out=folder/'scale_transfer';out.mkdir(exist_ok=True)
    with (out/'metrics.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    np.savez_compressed(out/'spectra.npz',time_seconds=[r['time_seconds'] for r in rows],k_upper=spectra[0]['k_upper'],**{k:np.stack([s[k] for s in spectra]) for k in ['energy_shell','advective_transfer_shell','advective_flux','stress_work_shell']})
    provenance={'method':'Even/odd free-slip channel reflection; Fourier 2/3 truncation; solenoidal projection of cell-centered density-weighted depth-mean velocity; dealiased resolved advective transfer and depth-dependent Reynolds-stress work separately.',
        'flux_sign':'Negative advective_flux means upscale; positive stress_work means energy entering the depth-mean flow from depth-dependent motions.',
        'large_scale_cutoff':'k < 0.75 * 2*pi/initial_seed_wavelength; excludes uniform mean; wavelengths greater than 1600 km for 1200-km seeds.',
        'units':'E: m^2/s^2 per unit reference column mass. T, flux, stress work: m^2/s^3. k: radians/m.',
        'limitations':'Diagnostic of projected, resolved mean flow; not an exact discrete WENO budget. Numerical dissipation, Coriolis transfer, sponge and pressure/discretization residuals are not closed here. Depth-dependent stress work includes all departures from the depth mean and is not uniquely buoyant convection. Use heated/control differences and budget/retained-energy metrics. A broad constant-flux inertial range is not guaranteed by this finite box.',
        'frame_count':len(rows),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'min_retained_energy_fraction':min(r['retained_energy_fraction'] for r in rows),'max_projection_removed_fraction':max(r['projection_removed_fraction_of_filtered_energy'] for r in rows),
        'status':'COMPLETE'}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Transfer diagnostics:',out,'frames:',len(rows),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path,nargs='?');p.add_argument('--self-test',action='store_true');a=p.parse_args()
    if a.self_test:self_test()
    if a.folder:analyze(a.folder)
    if not a.self_test and not a.folder:p.error('Provide a completed segment or --self-test')

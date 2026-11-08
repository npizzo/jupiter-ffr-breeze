"""Postprocess projected barotropic self-transfer and baroclinic-stress work.

Specific KE units are m²/s²; transfer is m²/s³. This is not a full 3-D or
exact discrete WENO budget. Negative cumulative self-flux denotes upscale transfer.
"""
import argparse,json,tomllib
from pathlib import Path
import h5py
import numpy as np


def extend(a, parity, channel):
    return np.concatenate((a,parity*a[::-1,:]),axis=0) if channel else a


def transfer(fields,lx,ly,channel=False):
    u=extend(fields['ubar'],1,channel);v=extend(fields['vbar'],-1,channel)
    ny,nx=u.shape;length_y=2*ly if channel else ly
    kx=2*np.pi*np.fft.fftfreq(nx,d=lx/nx)[None,:]
    ky=2*np.pi*np.fft.fftfreq(ny,d=length_y/ny)[:,None]
    k2=kx*kx+ky*ky;safe=np.where(k2==0,1,k2)
    mask=(abs(np.fft.fftfreq(nx)*nx)[None,:]<nx/3)&(abs(np.fft.fftfreq(ny)*ny)[:,None]<ny/3)
    fft=lambda a:np.fft.fft2(a,norm='forward')
    inv=lambda a:np.fft.ifft2(a,norm='forward').real
    rawu,rawv=fft(u),fft(v)
    uh,vh=rawu*mask,rawv*mask
    div=kx*uh+ky*vh
    pu=uh-kx*div/safe;pv=vh-ky*div/safe
    projection_fraction=float(np.sum(abs(pu-uh)**2+abs(pv-vh)**2)/max(np.sum(abs(uh)**2+abs(vh)**2),1e-300))
    u,v=inv(pu),inv(pv)
    nu=-fft(u*inv(1j*kx*pu)+v*inv(1j*ky*pu))*mask
    nv=-fft(u*inv(1j*kx*pv)+v*inv(1j*ky*pv))*mask
    work=np.real(np.conj(pu)*nu+np.conj(pv)*nv)
    scale=float(np.sum(abs(work)))
    residual=float(np.sum(work))
    rms=float(np.sqrt(np.mean(u*u+v*v)))
    roundoff_scale=rms**3*float(np.sqrt(k2).max())
    assert abs(residual)<1e-10*scale+1e-13*max(roundoff_scale,1e-30), 'Barotropic self-transfer does not conserve energy'
    rs=[fft(extend(fields[key],parity,channel))*mask for key,parity in [('Rxx',1),('Rxy',-1),('Ryy',1)]]
    ex=-1j*(kx*rs[0]+ky*rs[1]);ey=-1j*(kx*rs[1]+ky*rs[2])
    eddy=np.real(np.conj(pu)*ex+np.conj(pv)*ey)
    dk=min(2*np.pi/lx,2*np.pi/length_y)
    shell=np.rint(np.sqrt(k2)/dk).astype(int)
    bins=lambda a:np.bincount(shell.ravel(),weights=a.ravel())
    energy=bins((abs(pu)**2+abs(pv)**2)/2);self_transfer=bins(work)
    assert np.isclose(energy.sum(),np.mean(u*u+v*v)/2,rtol=1e-12,atol=1e-20)
    return {'k_rad_m':np.arange(len(energy))*dk,'energy':energy,'self_transfer':self_transfer,
            'self_flux':-np.cumsum(self_transfer),'baroclinic_stress_work':bins(eddy),
            'projection_energy_fraction':projection_fraction,'self_transfer_residual':residual,
            'removed_high_wavenumber_energy_fraction':float(np.sum((abs(rawu)**2+abs(rawv)**2)*(~mask))/max(np.sum(abs(rawu)**2+abs(rawv)**2),1e-300)),
            'zonal_mean_energy_fraction':float(np.mean(u.mean(axis=1)**2+v.mean(axis=1)**2)/max(np.mean(u*u+v*v),1e-300))}


def main():
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    m=tomllib.loads((args.case/'metadata.toml').read_text());channel=m['meridional_boundary']=='channel'
    records=[];times=[]
    for path in sorted(args.case.glob('frame_*.jld2')):
        with h5py.File(path) as f:
            fields={k:f[k][:] for k in ('ubar','vbar','Rxx','Rxy','Ryy')};times.append(float(f['t'][()]))
        assert all(np.isfinite(a).all() for a in fields.values())
        records.append(transfer(fields,m['lx'],m['ly'],channel))
    assert len(times)>1 and np.all(np.diff(times)>0)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output,time_seconds=times,**{key:np.asarray([r[key] for r in records]) for key in records[0]})
    args.output.with_suffix('.json').write_text(json.dumps({'case':m['name'],'frames':len(times),'method':'density-weighted full-depth horizontal velocity; 2/3 spectral truncation and horizontal solenoidal projection; even/odd extension for free-slip channel',
        'negative_self_flux':'net upscale transfer within the projected barotropic mode',
        'limits':'Excludes full 3-D transfer, pressure/buoyancy/viscous/numerical terms and exact WENO closure. Baroclinic-stress work is a separate diagnostic. Inspect projection/truncation losses and time averages before interpretation.'},indent=2)+'\n')
    print('Saved barotropic transfer diagnostics:',args.output)

if __name__=='__main__':main()

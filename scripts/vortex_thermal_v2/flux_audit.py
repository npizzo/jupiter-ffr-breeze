"""Additional read-only separation of vortex mean transport and convective covariance."""
import sys,csv,json
from pathlib import Path
import numpy as np,h5py
from analyze import core,low_vorticity,periodic_components,writecsv
CP=2.5
out=Path(sys.argv[1]);assert (out/'ANALYSIS_PASS.json').is_file()
s=json.loads((out/'summary.json').read_text());result={};allrows=[]
for case,prov in s['provenance'].items():
    fields={};frames={}
    for seg in map(Path,prov['segments']):
        for p in (seg/'fields3d').glob('*.jld2'):
            with h5py.File(p) as h:t=float(h['t'][()])
            if t==0 or t>=380:fields[t]=p
        for p in seg.glob('frame_*.jld2'):
            with h5py.File(p) as h:t=float(h['t'][()])
            if t in fields:frames[t]=p
    tracks={(float(r['time']),r['region']):r for r in csv.DictReader((out/f'{case}_vortex_tracks.csv').open())}
    init={};rows=[];maxerr=0.
    for t,p in sorted(fields.items()):
        if t not in frames:continue
        with h5py.File(frames[t]) as h:ub=h['ubar'][()].astype(float);vb=h['vbar'][()].astype(float)
        zb,zl,kr=low_vorticity(ub,vb);vort=zb if t==0 else zl;regions={}
        for sign,region in [(1,'cyclone'),(-1,'anticyclone')]:
            if t==0:m,rec=core(vort,sign)
            else:
                row=tracks[t,region];labs=periodic_components(sign*vort>.4*np.max(sign*vort));candidates=[]
                for lab in np.unique(labs):
                    if not lab:continue
                    mask=labs==lab;n=mask.shape[0];coords=(np.arange(n)+.5)*4/n
                    def center(weights):return float(np.mod(np.angle(np.sum(weights*np.exp(2j*np.pi*coords/4))),2*np.pi)*4/(2*np.pi))
                    x,y=center(mask.sum(0)),center(mask.sum(1));dx=abs(x-float(row['x']));dy=abs(y-float(row['y']));dist=min(dx,4-dx)**2+min(dy,4-dy)**2
                    candidates.append((dist,mask))
                dist,m=min(candidates,key=lambda q:q[0]);assert dist<1e-10,(case,t,region,dist)
            regions[region]=m
        regions['background']=~(regions['cyclone']|regions['anticyclone'])
        with h5py.File(p) as h:rho=h['density'][()].astype(float);T=h['temperature'][()].astype(float);wf=h['w_face'][()].astype(float);z=h['z'][()]
        w=.5*(wf[:-1]+wf[1:]);hp=lambda a:np.fft.ifft2(np.fft.fft2(a,axes=(-2,-1))*(kr>2),axes=(-2,-1)).real
        wh=hp(w);Th=hp(T);Ta=T-T.mean((1,2))[:,None,None]
        for reg,m in regions.items():
            rp=rho[:,m];wp=w[:,m];Tp=T[:,m];rm=rp.mean(1);wt=(rp*wp).mean(1)/rm;tt=(rp*Tp).mean(1)/rm
            raw=CP*(rp*wp*Tp).mean(1);meanterm=CP*rm*wt*tt;cov=CP*(rp*(wp-wt[:,None])*(Tp-tt[:,None])).mean(1)
            err=float(np.max(abs(raw-meanterm-cov)));assert err<1e-11;maxerr=max(maxerr,err)
            hpf=CP*(rp*wh[:,m]*Th[:,m]).mean(1)
            for k,zz in enumerate(z):rows.append(dict(case=case,time=t,region=reg,z=float(zz),temperature_anomaly=float(Ta[:,m].mean(1)[k]),raw_enthalpy_flux=float(raw[k]),conditional_mean_transport=float(meanterm[k]),within_region_covariance_flux=float(cov[k]),smallscale_enthalpy_flux=float(hpf[k])))
            if t==0:init[reg]={'temperature_anomaly_bulk':float(Ta[:,m].mean(1)[z<.78].mean()),'w_rms':float(np.sqrt(np.mean(w[:,m]**2)))}
        print('FLUX_AUDIT',case,t,flush=True)
    def avg(a,b):
        return {reg:{key:float(np.mean([r[key] for r in rows if a<=r['time']<=b and r['region']==reg and r['z']<.78])) for key in ['within_region_covariance_flux','smallscale_enthalpy_flux','conditional_mean_transport']} for reg in regions}
    last=max(fields);result[case]={'initial_full_seed_vorticity_conditioned':init,'matched_400_480':avg(400,480),'last_100':avg(last-100,last),'decomposition_error_max':maxerr};allrows+=rows
writecsv(out/'flux_decomposition.csv',allrows)
(out/'flux_audit.json').write_text(json.dumps({'status':'PASS','cases':result,'definitions':{'within_region_covariance_flux':'Cp mean_region[rho(w-Favre_region_mean(w))(T-Favre_region_mean(T))] at each height','smallscale_enthalpy_flux':'Cp mean_region[rho w_highpass T_highpass], horizontal Fourier |k|>2; cutoff-dependent diagnostic','initial':'At t=0 use full depth-mean vorticity of seeded eddies; later use tracked low-pass LSV cores. These are different spatial masks, not a Lagrangian temperature comparison.'},'limitation':'Conditional associations and scale-filtered fluxes do not close a vortex energy budget or demonstrate causality.'},indent=2)+'\n')
print('FLUX_AUDIT_COMPLETE',json.dumps(result),flush=True)

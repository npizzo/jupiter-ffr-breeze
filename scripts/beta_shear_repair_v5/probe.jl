function add_probe!(m,p,cfg)
    amplitude=get(cfg,"probe_amplitude",0.)
    amplitude==0 && return
    a=arrays(m)
    # Smooth mass-flux streamfunction on y/z corners; exactly zero normal wall flow.
    lo=get(cfg,"probe_low",0.);hi=get(cfg,"probe_high",.65)
    psi=[exp(-((y-0.)/.75)^2)*sinpi((y+p.ly/2)/p.ly)^4*(lo<z<hi ? sinpi((z-lo)/(hi-lo))^4 : 0.)*cos(2pi*y/.75) for y=range(-p.ly/2,p.ly/2;length=p.ny+1),z=range(0,1;length=p.nz+1)]
    psi[[1,end],:].=0.;psi[:,[1,end]].=0.
    mv=diff(psi;dims=2)./p.dz;mw=-diff(psi;dims=1)./p.dy
    ry=yface_density(a.rho)[1,:,:]
    rz=cat(a.rho[1,:,1:1],.5 .*(a.rho[1,:,1:end-1].+a.rho[1,:,2:end]),a.rho[1,:,end:end];dims=2)
    dv=mv./ry;dw=mw./rz
    norm=sqrt((sum(dv.^2)+sum(dw.^2))/(p.ny*p.nz))
    dv.*=amplitude/norm;dw.*=amplitude/norm
    vv=a.vf.+reshape(dv,1,p.ny+1,p.nz);ww=a.wf.+reshape(dw,1,p.ny,p.nz+1)
    set!(m;u=a.uf,v=vv,w=ww)
    Oceananigans.TimeSteppers.update_state!(m)
    @assert maximum(abs,Array(interior(m.dynamics.dry_density)).-a.rho)<1e-12
    @assert maximum(abs,Array(interior(m.temperature)).-a.T)<1e-12
    divergence=diff(mv;dims=1)./p.dy.+diff(mw;dims=2)./p.dz
    @assert maximum(abs,divergence)<1e-10
    @assert maximum(abs,vv[:,[1,end],:])==0 && maximum(abs,ww[:,:,[1,end]])==0
    p.seed_info["probe_rms"]=amplitude
    p.seed_info["probe_description"]="Smooth x-uniform divergence-free mass-flux streamfunction, fixed physical wavelengths; y/z face velocity RMS normalization."
    nothing
end

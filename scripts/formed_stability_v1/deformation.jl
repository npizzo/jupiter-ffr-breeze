using Oceananigans.Solvers: FFTW
function deform!(m,rho,cfg,out)
    n,ny,nz=cfg["nx"],cfg["ny"],cfg["nz"];Lx,Ly=cfg["lx"],cfg["ly"];dx,dy=Lx/n,Ly/ny
    h=reshape(rho,1,1,:);a=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(Oceananigans.prognostic_fields(m)))
    u=a["ρu"] ./ h;v=a["ρv"] ./ h;ub=dropdims(sum(h .* u;dims=3) ./ sum(h);dims=3);vb=dropdims(sum(h .* v;dims=3) ./ sum(h);dims=3)
    if haskey(cfg,"source_state")
        old=load(cfg["source_state"])
        for (key,q) in [("u_face",u),("v_face",v),("w_face",Array(interior(m.velocities.w))),("dye_density",a["dye"])]
            @assert maximum(abs,q .- old[key])<1e-10max(maximum(abs,q),1.) "Actual formed-vortex checkpoint differs from saved fields"
        end
        @assert isapprox(m.clock.time,100*86400.;atol=1e-8)
    end
    uc,vc=nb_center(ub,vb);zz=nb_ccurl(uc,vc,dx,dy)
    # Smoothing by a physical Gaussian multiplier with an even channel extension.
    ext=hcat(zz,reverse(zz;dims=2));nx,nn=size(ext)
    kx=reshape([min(i,nx-i)*2pi/Lx for i=0:nx-1],:,1);ky=reshape([min(j,nn-j)*pi/Ly for j=0:nn-1],1,:)
    filtered=real.(FFTW.ifft(FFTW.fft(ext) .* exp.(-.5*(.08min(Lx,Ly))^2 .* (kx .^ 2  .+ ky .^ 2))))[:,1:ny]
    ix,jy=Tuple(argmax(filtered));cx=(ix-.5)*dx;cy=(jy-.5)*dy-Ly/2
    x=reshape((0:n-1) .* dx,:,1);y=reshape((0:ny) .* dy .- Ly/2,1,:)
    X=Lx/pi .* sinpi.((x .- cx) ./ Lx);Y=y .- cy;radius=.15min(Lx,Ly)
    psi=(X .^ 2  .- Y .^ 2) ./ radius^2 .* exp.(-.5 .* (X .^ 2  .+ Y .^ 2) ./ radius^2) .* sinpi.((y .+ Ly/2) ./ Ly) .^ 4
    psi[:,[1,end]].=0
    du=-diff(psi;dims=2) ./ dy;dv=(circshift(psi,(-1,0)) .- psi) ./ dx
    speed=sqrt(mean(ub .^ 2)+mean(.5 .* (vb[:,1:end-1] .^ 2  .+ vb[:,2:end] .^ 2)))
    perturb_rms=sqrt(mean(du .^ 2)+mean(.5 .* (dv[:,1:end-1] .^ 2  .+ dv[:,2:end] .^ 2)))
    eps=cfg["deformation_fraction"];factor=eps*speed/perturb_rms;du .*=factor;dv .*=factor;psi .*=factor
    @assert maximum(abs,(circshift(du,(-1,0)) .- du) ./ dx .+ diff(dv;dims=2) ./ dy)<1e-12
    @assert maximum(abs,dv[:,[1,end]])==0
    if eps!=0
        set!(m.momentum.ρu,a["ρu"] .+ h .* reshape(du,n,ny,1));set!(m.momentum.ρv,a["ρv"] .+ h .* reshape(dv,n,ny+1,1))
        pc=.25 .* (psi[:,1:end-1] .+ psi[:,2:end] .+ circshift(psi,(-1,0))[:,1:end-1] .+ circshift(psi,(-1,0))[:,2:end])
        f=2cfg["rotation_rate"]*sind(cfg["latitude"])
        set!(m.dynamics.pressure_anomaly,Array(interior(m.dynamics.pressure_anomaly)) .+ reshape(f .* pc,n,ny,1))
    end
    b=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(Oceananigans.prognostic_fields(m)))
    for k in keys(a);k in ["ρu","ρv"] && continue;@assert a[k]==b[k];end
    div=mass_divergence(b,dx,dy,cfg["height"]/nz);@assert maximum(abs,div)<1e-10
    jldsave(joinpath(out,"deformation.jld2");du,dv,psi,center=[cx,cy],radius,fraction=eps,parent_rms=speed,parent_vorticity=zz,filtered_parent=filtered)
    report=Dict("status"=>"PASS","fraction"=>eps,"parent_barotropic_rms"=>speed,"core_center"=>[cx,cy],"radius"=>radius,"max_mass_divergence"=>maximum(abs,div),"vertical_velocity_and_thermodynamics_unchanged"=>true,"pressure"=>"Geostrophic increment; projection resolves remaining nonlinear adjustment; analyze after first inertial period")
    open(joinpath(out,"PERTURBATION_PASS.toml"),"w") do io;TOML.print(io,report);end
end

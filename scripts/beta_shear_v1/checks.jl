using Oceananigans.BoundaryConditions: fill_halo_regions!
function operator_checks(cfg,out)
    c=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>24,"initial_rms"=>0.,"jet_speed"=>0.))
    m,p=build_model(c;arch=CPU(),perturb=false)
    T=[reference(zcell(k,p)).T+.03cospi((j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    Q=[conduction(i,j,k,T,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    qt=mean(2p.kc .* (T[:,:,end].-p.tt)./p.dz)
    fluxerr=abs(mean(Q)-(p.flux-qt));@assert fluxerr<1e-11
    rho=[reference(zcell(k,p)).rho for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(p.nu,(x,y,z)->.001*(1+.2cospi(y/4)));fill_halo_regions!(p.nu)
    F=(T=T,ρᵈ=rho)
    qs=[sgs_heating(i,j,k,F,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    sgs_err=abs(mean(qs));@assert sgs_err<1e-10
    M=(ρu=ones(24,25,32),ρv=2ones(24,25,32),ρw=3ones(24,25,32))
    fy=2p.omega*cosd(c["latitude"]);err=0.
    for j=3:p.ny-2
        yc=(j-.5)*p.dy;yf=(j-1)*p.dy
        err=max(err,abs(Oceananigans.Coriolis.x_f_cross_U(6,j,6,m.grid,m.coriolis,M)-(3fy-2*(p.f+p.beta*(yc-p.ly/2)))))
        err=max(err,abs(Oceananigans.Coriolis.y_f_cross_U(6,j,6,m.grid,m.coriolis,M)-(p.f+p.beta*(yf-p.ly/2))))
        err=max(err,abs(Oceananigans.Coriolis.z_f_cross_U(6,j,6,m.grid,m.coriolis,M)+fy))
    end
    @assert err<1e-12
    # A unit northward current must give vorticity tendency -beta.
    Mv=(ρu=zeros(24,25,32),ρv=ones(24,25,32),ρw=zeros(24,25,32))
    cx(j)=Oceananigans.Coriolis.x_f_cross_U(6,j,6,m.grid,m.coriolis,Mv)
    betacurl=(cx(7)-cx(6))/p.dy
    @assert abs(betacurl+p.beta)<1e-12
    write_toml(joinpath(out,"OPERATORS_"*cfg["name"]*".toml"),Dict("status"=>"PASS","closed_channel_conduction_error"=>fluxerr,"closed_channel_sgs_error"=>sgs_err,"coriolis_error"=>err,"uniform_v_curl_tendency"=>betacurl,"expected_curl"=>-p.beta))
end
function balance_check(m,p)
    a=arrays(m);P=Array(interior(m.dynamics.pressure));Pr=Array(interior(m.dynamics.reference_state.pressure));r0=Array(interior(m.dynamics.reference_state.density))
    mx=Array(interior(m.momentum.ρu));rx=.5 .*(a.rho.+circshift(a.rho,(1,0,0)))
    iz=.25 .*(mx[:,:,1:end-1].+mx[:,:,2:end].+circshift(mx,(-1,0,0))[:,:,1:end-1].+circshift(mx,(-1,0,0))[:,:,2:end])
    rz=.5 .*(a.rho[:,:,1:end-1].+a.rho[:,:,2:end]);fy=2p.omega*cosd(67.5)
    az=(-diff(P.-Pr;dims=3)./p.dz.-GRAV.*(rz.-.5 .*(r0[:,:,1:end-1].+r0[:,:,2:end])).+fy.*iz)./rz
    scale=abs(p.f)*sqrt(mean(a.u.^2+a.v.^2));err=maximum(abs,az)/scale;@assert err<1e-8
    # Report horizontal acceleration mismatch in the interior; beta is not hidden by rebalancing IC.
    g=m.grid;M=m.momentum;res=0.;count=0
    uc=a.u;vc=a.v
    advx=uc.*(circshift(uc,(-1,0,0)).-circshift(uc,(1,0,0)))./(2p.dx).+vc.*(circshift(uc,(0,-1,0)).-circshift(uc,(0,1,0)))./(2p.dy)
    advy=uc.*(circshift(vc,(-1,0,0)).-circshift(vc,(1,0,0)))./(2p.dx).+vc.*(circshift(vc,(0,-1,0)).-circshift(vc,(0,1,0)))./(2p.dy)
    ry=yface_density(a.rho)
    ax=zeros(size(a.uf));ay=zeros(size(a.vf))
    for k=2:p.nz-1,j=2:p.ny-1,i=1:p.nx
        im=i==1 ? p.nx : i-1
        ax[i,j,k]=-(P[i,j,k]-P[im,j,k])/(p.dx*rx[i,j,k])-Oceananigans.Coriolis.x_f_cross_U(i,j,k,g,m.coriolis,M)/rx[i,j,k]
        ay[i,j,k]=-(P[i,j,k]-P[i,j-1,k])/(p.dy*ry[i,j,k])-Oceananigans.Coriolis.y_f_cross_U(i,j,k,g,m.coriolis,M)/ry[i,j,k]
    end
    acx=.5 .*(ax.+circshift(ax,(-1,0,0)));acy=.5 .*(ay[:,1:end-1,:].+ay[:,2:end,:]);sel=3:p.ny-2
    h=sqrt(mean((acx[:,sel,2:end-1].-advx[:,sel,2:end-1]).^2+(acy[:,sel,2:end-1].-advy[:,sel,2:end-1]).^2))/scale
    @assert h<.2 "Channel seed imbalance must be inspected before production"
    if get(p.seed_info,"initial_rms",0.)<1e-3;@assert h<.02 "Jet-only balance must be tighter than the turbulent seed";end
    mach=maximum(sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T))
    @assert mach<.5 "Initial jet-plus-eddy state must remain subsonic below Mach 0.5"
    jmid=div(p.ny,2)
    shear=(jet_profile(p.dy/2,p)-jet_profile(-p.dy/2,p))/p.dy
    @assert p.jet_speed==0 || shear>0
    @assert jet_profile(-p.ly/2,p)==0 && jet_profile(p.ly/2,p)==0
    Dict("initial_max_mach"=>mach,"central_jet_du_dy"=>shear,"vertical_imbalance_over_fU"=>err,"interior_horizontal_imbalance_over_fU"=>h,"normal_wall_velocity"=>maximum(abs,a.vf[:,[1,end],:]))
end

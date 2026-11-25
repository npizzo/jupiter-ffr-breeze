using LinearAlgebra
# A common monochromatic streamfunction gives mixed signs without favoring a core.
# All four integer modes have |n|²=5 and change sign under a (L/2,L/2) shift.
function seed_velocity(p,cfg)
    phase=[0.31,1.17,2.03,2.71];modes=[(1,2),(2,1),(-1,2),(-2,1)]
    psi=zeros(p.nx,p.ny)
    for ((kx,ky),ph) in zip(modes,phase),j=1:p.ny,i=1:p.nx
        psi[i,j]+=cos(2pi*(kx*(i-1)/p.nx+ky*(j-1)/p.ny)+ph)
    end
    psi.-=mean(psi)
    u=-(circshift(psi,(0,-1)).-psi)./p.dy
    v=(circshift(psi,(-1,0)).-psi)./p.dx
    fac=cfg["initial_rms"]/sqrt(mean(u.^2+v.^2));psi.*=fac;u.*=fac;v.*=fac
    return psi,u,v
end
function density_to_pressure(rho,Pref,rhoref,u,anchor,p,latitude)
    # Integrate the native vertical balance instead of inverting face averages.
    # Inverting face averages admits a grid-scale alternating density mode.
    fy=2p.omega*cosd(latitude)
    rx=.5 .* (rho.+circshift(rho,(1,0,0)))
    mx=rx.*reshape(u,p.nx,p.ny,1)
    rz=.5 .* (rho[:,:,1:end-1].+rho[:,:,2:end])
    iz=.25 .* (mx[:,:,1:end-1].+mx[:,:,2:end].+circshift(mx,(-1,0,0))[:,:,1:end-1].+circshift(mx,(-1,0,0))[:,:,2:end])
    rref=reshape(.5 .* (rhoref[1:end-1].+rhoref[2:end]),1,1,:)
    difference=diff(Pref;dims=3).+p.dz.*(-GRAV.*(rz.-rref).+fy.*iz)
    P=zeros(size(rho));P[:,:,end].=anchor
    for k=p.nz-1:-1:1;P[:,:,k].=P[:,:,k+1].-difference[:,:,k];end
    P
end
function seed_moments(rho,u,v,p)
    rx=.5 .* (rho.+circshift(rho,(1,0,0)));ry=.5 .* (rho.+circshift(rho,(0,1,0)))
    dv=p.dx*p.dy*p.dz
    [sum(rho)*dv, .5sum(rx.*reshape(u,p.nx,p.ny,1).^2+ry.*reshape(v,p.nx,p.ny,1).^2)*dv]
end
function apply_seed!(m,p,cfg;perturb=true)
    get(cfg,"initial_rms",0.)>0||return Dict{String,Any}()
    psi,u,v=seed_velocity(p,cfg);z=[zcell(k,p) for k=1:p.nz]
    Pref=repeat(reshape(vec(Array(interior(m.dynamics.reference_state.pressure))[:, :, :][1,1,:]),1,1,:),p.nx,p.ny,1)
    rhoref=vec(Array(interior(m.dynamics.reference_state.density))[:, :, :][1,1,:])
    Rref=repeat(reshape(rhoref,1,1,:),p.nx,p.ny,1)
    pc=.25 .* (psi.+circshift(psi,(-1,0)).+circshift(psi,(0,-1)).+circshift(psi,(-1,-1)))
    uc=.5 .* (u.+circshift(u,(-1,0)));vc=.5 .* (v.+circshift(v,(0,-1)))
    kh2=(2pi/4)^2*5
    Q=p.f.*pc.-.5 .*(uc.^2 .+vc.^2 .+kh2 .*pc.^2);Q.-=mean(Q)
    P=[Pref[i,j,k]*reference(z[k]-Q[i,j]/GRAV).p/reference(z[k]).p for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    # Smooth density from the displaced polytrope and the tilted rotation balance.
    factor=GRAV./(GRAV.-2p.omega*cosd(cfg["latitude"]).*uc)
    rho=[Rref[i,j,k]*reference(z[k]-Q[i,j]/GRAV).rho/reference(z[k]).rho*factor[i,j] for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    target=seed_moments(Rref,u,v,p)
    ke_shape=.5 .*(uc.^2 .+vc.^2);ke_shape=(ke_shape.-mean(ke_shape))./cfg["initial_rms"]^2
    r1=Rref;r2=Rref.*reshape(ke_shape,p.nx,p.ny,1)
    A=hcat(seed_moments(r1,u,v,p),seed_moments(r2,u,v,p));b=target.-seed_moments(rho,u,v,p)
    coefficients=(A./reshape(target,2,1))\(b./target)
    @assert maximum(abs,coefficients)<.05 "Balance mass/KE correction unexpectedly large"
    rho.+=coefficients[1].*r1.+coefficients[2].*r2
    anchor=P[:,:,end].+Pref[:,:,end].*(coefficients[1].+coefficients[2].*ke_shape)
    P=density_to_pressure(rho,Pref,rhoref,u,anchor,p,cfg["latitude"])
    zz=reshape(z,1,1,:);dv=p.dx*p.dy*p.dz
    Eref=sum(CV.*Pref.+GRAV.*Rref.*zz)*dv
    # Rotation-dependent balanced potential/internal energy is reported, not removed
    # by changing the entropy profile to enforce an artificial common total energy.
    pressure_offset=0.
    @assert minimum(P)>0 && minimum(rho)>0
    @assert maximum(abs,rho./Rref.-1)<.75 "Seed is too thermodynamically disruptive"
    rng=MersenneTwister(cfg["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
    # The same zero-horizontal-mean absolute pressure noise adds zero integrated heat.
    if perturb
        P.+=(5/3)*cfg["perturbation"].*Pref.*noise.*reshape(sinpi.(clamp.(z./DC,0.,1.)),1,1,:)
    end
    theta=P.^.6 ./rho
    set!(m;θ=theta,ρ=rho,qᵗ=0.,u=repeat(reshape(u,p.nx,p.ny,1),1,1,p.nz),v=repeat(reshape(v,p.nx,p.ny,1),1,1,p.nz),w=0.)
    Oceananigans.TimeSteppers.update_state!(m)
    actual=arrays(m);got=seed_moments(actual.rho,u,v,p)
    @assert maximum(abs,got./target.-1)<1e-10
    zeta=(v.-circshift(v,(1,0)))./p.dx.-(u.-circshift(u,(0,1)))./p.dy
    div=(circshift(u,(-1,0)).-u)./p.dx.+(circshift(v,(0,-1)).-v)./p.dy
    @assert maximum(abs,div)<1e-12
    @assert abs(mean(zeta))<1e-12
    @assert abs(mean(zeta.^3)/mean(zeta.^2)^1.5)<1e-10
    pressure_error=maximum(abs,Array(interior(m.dynamics.pressure)).-P)/maximum(P)
    @assert pressure_error<1e-10
    Etarget=sum(CV.*P.+GRAV.*rho.*zz)*dv+target[2]
    Eactual=sum(actual.rho.*(CV.*actual.T.+GRAV.*zz))*dv+got[2]
    @assert abs(Eactual/Etarget-1)<1e-10
    if cfg["formulation"]=="TotalEnergy"
        Enative=sum(interior(m.formulation.energy_density))*dv
        @assert abs(Enative/Eactual-1)<1e-10 "Seed total energy does not match physical state"
    end
    @assert maximum(abs,diff(actual.uf;dims=3))<1e-12
    @assert maximum(abs,diff(actual.vf;dims=3))<1e-12
    @assert maximum(abs,actual.wf)<1e-12

    N2=GRAV.*diff(log.(P.^.6 ./rho);dims=3)./p.dz
    zface=.5 .* (z[1:end-1].+z[2:end]);deep=findall(zface.<.65);cap=findall(zface.>.85)
    stability=Dict("N2_deep_min"=>minimum(N2[:,:,deep]),"N2_deep_max"=>maximum(N2[:,:,deep]),"N2_cap_min"=>minimum(N2[:,:,cap]),"N2_cap_max"=>maximum(N2[:,:,cap]),"balanced_energy_difference"=>Eactual-(Eref+target[2]))
    println("INITIAL_STABILITY ",stability);flush(stdout)
    @assert stability["N2_cap_min"]>0 "Balanced eddies must preserve a stable upper cap"
    @assert max(abs(stability["N2_deep_min"]),abs(stability["N2_deep_max"]))<.2 "Seed spuriously changes bulk convective stability"
    initial_info=Dict("status"=>"INITIALIZED","initial_rms"=>cfg["initial_rms"],"wavelength"=>4/sqrt(5),"turnover_time"=>4/(sqrt(5)*cfg["initial_rms"]),"mass"=>got[1],"kinetic_energy"=>got[2],"total_energy"=>Eactual,"native_vorticity_skewness"=>mean(zeta.^3)/mean(zeta.^2)^1.5,"max_abs_zeta_over_f"=>maximum(abs,zeta)/abs(p.f),"velocity_divergence_max"=>maximum(abs,div),"pressure_readback_error"=>pressure_error,"density_relative_departure_max"=>maximum(abs,rho./Rref.-1),"matching_coefficients"=>collect(coefficients),"pressure_offset"=>pressure_offset,"thermal_equilibrium"=>false,"note"=>"Finite-amplitude mixed-sign flow; mass and kinetic energy matched; pressure/density and balanced thermal/potential energy differ for rotation-dependent balance. No mature thermal parent.")
    merge!(initial_info,stability)
    return initial_info
end
function initial_balance(m,p)
    a=arrays(m);u=a.uf;v=a.vf;rho=a.rho
    P=Array(interior(m.dynamics.pressure));Pref=Array(interior(m.dynamics.reference_state.pressure));r0=Array(interior(m.dynamics.reference_state.density))
    rx=.5 .*(rho.+circshift(rho,(1,0,0)));ry=.5 .*(rho.+circshift(rho,(0,1,0)))
    mx=Array(interior(m.momentum.ρu));my=Array(interior(m.momentum.ρv))
    fz=p.f;fy=2p.omega*cosd(67.5)
    iv=.25 .*(my.+circshift(my,(1,0,0)).+circshift(my,(0,-1,0)).+circshift(my,(1,-1,0)))
    iu=.25 .*(mx.+circshift(mx,(-1,0,0)).+circshift(mx,(0,1,0)).+circshift(mx,(-1,1,0)))
    # Horizontal pressure+Coriolis accelerations; compare against centered nonlinear advection.
    ax=-(P.-circshift(P,(1,0,0)))./(p.dx.*rx).+fz.*iv./rx
    ay=-(P.-circshift(P,(0,1,0)))./(p.dy.*ry).-fz.*iu./ry
    uc=a.u;vc=a.v
    advx=uc.*(circshift(uc,(-1,0,0)).-circshift(uc,(1,0,0)))./(2p.dx).+vc.*(circshift(uc,(0,-1,0)).-circshift(uc,(0,1,0)))./(2p.dy)
    advy=uc.*(circshift(vc,(-1,0,0)).-circshift(vc,(1,0,0)))./(2p.dx).+vc.*(circshift(vc,(0,-1,0)).-circshift(vc,(0,1,0)))./(2p.dy)
    acx=.5 .*(ax.+circshift(ax,(-1,0,0)));acy=.5 .*(ay.+circshift(ay,(0,-1,0)))
    scale=abs(fz)*sqrt(mean(uc.^2 .+vc.^2))
    horizontal_rms=sqrt(mean((acx.-advx).^2 .+(acy.-advy).^2))/scale
    iz=.25 .*(mx[:,:,1:end-1].+mx[:,:,2:end].+circshift(mx,(-1,0,0))[:,:,1:end-1].+circshift(mx,(-1,0,0))[:,:,2:end])
    rz=.5 .*(rho[:,:,1:end-1].+rho[:,:,2:end])
    az=(-(diff(P.-Pref;dims=3))./p.dz.-GRAV.*(rz.-.5 .*(r0[:,:,1:end-1].+r0[:,:,2:end])).+fy.*iz)./rz
    return Dict("horizontal_imbalance_rms_over_fU"=>horizontal_rms,"vertical_imbalance_max_over_fU"=>maximum(abs,az)/scale,"coriolis_work_over_bottom_power"=>sum(fz.*u.*iv.-fz.*v.*iu)*p.dx*p.dy*p.dz/(16p.flux))
end

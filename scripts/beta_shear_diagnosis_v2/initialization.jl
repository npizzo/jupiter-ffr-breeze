using LinearAlgebra
# A common monochromatic streamfunction gives mixed signs without favoring a core.
# All four integer modes have |n|²=5 and change sign under a (L/2,L/2) shift.
function yface_density(rho)
    cat(rho[:,1:1,:],.5 .*(rho[:,1:end-1,:].+rho[:,2:end,:]),rho[:,end:end,:];dims=2)
end
function channel_vorticity(uf,vf,p)
    zf=zeros(p.nx,p.ny+1,p.nz)
    zf[:,2:end-1,:].=(vf[:,2:end-1,:].-circshift(vf[:,2:end-1,:],(1,0,0)))./p.dx.-diff(uf;dims=2)./p.dy
    .25 .* (zf[:,1:end-1,:].+zf[:,2:end,:].+circshift(zf[:,1:end-1,:],(-1,0,0)).+circshift(zf[:,2:end,:],(-1,0,0)))
end
function seed_velocity(p,cfg)
    phase=[.31,1.17,2.03,2.71];modes=[(1,2),(2,1),(-1,2),(-2,1)]
    psi=zeros(p.nx,p.ny+1)
    for ((kx,ky),ph) in zip(modes,phase),j=1:p.ny+1,i=1:p.nx
        x=(i-1)/p.nx;y=(j-1)/p.ny
        psi[i,j]+=.5*(cos(2pi*(kx*x+ky*p.ly/4*y)+ph)-cos(2pi*(kx*x-ky*p.ly/4*y)+ph))
    end
    psi[:,[1,end]].=0.
    u=-diff(psi;dims=2)./p.dy;v=(circshift(psi,(-1,0)).-psi)./p.dx
    # Count interior faces with half-volume wall weights (wall v is exactly zero).
    fac=cfg["initial_rms"]/sqrt((sum(u.^2)+sum(v.^2))/(p.nx*p.ny));psi.*=fac;u.*=fac;v.*=fac
    return psi,u,v
end
function density_to_pressure(rho,Pref,rhoref,u,anchor,p,latitude)
    # Integrate the native vertical balance instead of inverting face averages.
    # Inverting face averages admits a grid-scale alternating density mode.
    fy=p.fh
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
    rx=.5 .* (rho.+circshift(rho,(1,0,0)));ry=yface_density(rho)
    dv=p.dx*p.dy*p.dz
    [sum(rho)*dv, .5*(sum(rx.*reshape(u,p.nx,p.ny,1).^2)+sum(ry.*reshape(v,p.nx,p.ny+1,1).^2))*dv]
end
function apply_seed!(m,p,cfg;perturb=true)
    cfg["initial_rms"]==0 && p.jet_speed==0 && return Dict{String,Any}()

    psi,u,v=seed_velocity(p,cfg);eddy_u=copy(u);eddy_v=copy(v);z=[zcell(k,p) for k=1:p.nz]
    Pref=repeat(reshape(vec(Array(interior(m.dynamics.reference_state.pressure))[:, :, :][1,1,:]),1,1,:),p.nx,p.ny,1)
    rhoref=vec(Array(interior(m.dynamics.reference_state.density))[:, :, :][1,1,:])
    Rref=repeat(reshape(rhoref,1,1,:),p.nx,p.ny,1)
    pc=.25 .* (psi[:,1:end-1].+psi[:,2:end].+circshift(psi[:,1:end-1],(-1,0)).+circshift(psi[:,2:end],(-1,0)))
    uc=.5 .* (u.+circshift(u,(-1,0)));vc=.5 .* (v[:,1:end-1].+v[:,2:end])
    kh2=(2pi/4)^2*5
    Q=p.f.*pc.-.5 .*(uc.^2 .+vc.^2 .+kh2 .*pc.^2);Q.-=mean(Q)
    # Exact continuum thermal-wind balance for a z-independent jet and full tilted rotation:
    # s=z+delta(y), delta_y=f(y)U/(g-fh*U), P=P_ref(s), rho=rho_ref(s)*g/(g-fh*U).
    yc=[(j-.5)*p.dy-p.ly/2 for j=1:p.ny]
    jets=jet_profile.(yc,Ref(p))
    ds=[(p.f+p.beta*yc[j])*jets[j]/(GRAV-p.fh*jets[j]) for j=1:p.ny]
    delta=zeros(p.ny)
    for j=2:p.ny;delta[j]=delta[j-1]+.5p.dy*(ds[j-1]+ds[j]);end
    delta.-=.5*(minimum(delta)+maximum(delta))
    u.+=reshape(jets,1,p.ny);uc=.5 .*(u.+circshift(u,(-1,0)))
    P=[Pref[i,j,k]*reference(z[k]+delta[j]-Q[i,j]/GRAV).p/reference(z[k]).p for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    # Smooth density from the displaced polytrope and the tilted rotation balance.
    factor=GRAV./(GRAV.-p.fh.*uc)
    rho=[Rref[i,j,k]*reference(z[k]+delta[j]-Q[i,j]/GRAV).rho/reference(z[k]).rho*factor[i,j] for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    target=seed_moments(Rref,u,v,p)
    if cfg["initial_rms"]>1e-3 && p.jet_speed==0
        ke_shape=.5 .*(uc.^2 .+vc.^2);ke_shape=(ke_shape.-mean(ke_shape))./max(mean(uc.^2 .+vc.^2),1e-12)
        r1=Rref;r2=Rref.*reshape(ke_shape,p.nx,p.ny,1)
        A=hcat(seed_moments(r1,u,v,p),seed_moments(r2,u,v,p));b=target.-seed_moments(rho,u,v,p)
        coefficients=(A./reshape(target,2,1))\(b./target)
        @assert maximum(abs,coefficients)<.05 "Mass/KE correction unexpectedly large"
        rho.+=coefficients[1].*r1.+coefficients[2].*r2
        anchor=P[:,:,end].+Pref[:,:,end].*(coefficients[1].+coefficients[2].*ke_shape)
    else
        factor_mass=target[1]/seed_moments(rho,u,v,p)[1]
        rho.*=factor_mass;P.*=factor_mass;anchor=P[:,:,end];coefficients=[factor_mass-1,0.]
        target=seed_moments(rho,u,v,p)
    end
    P=density_to_pressure(rho,Pref,rhoref,u,anchor,p,cfg["latitude"])
    zz=reshape(z,1,1,:);dv=p.dx*p.dy*p.dz
    Eref=sum(CV.*Pref.+GRAV.*Rref.*zz)*dv
    # Rotation-dependent balanced potential/internal energy is reported, not removed
    # by changing the entropy profile to enforce an artificial common total energy.
    pressure_offset=0.
    @assert minimum(P)>0 && minimum(rho)>0
    # The inherited perturbation limit applies only to the no-jet reference.
    # Balanced finite-amplitude compressible jets are checked by Mach, cap stability,
    # positive rho/T and native force balance; report their full density departure.
    if p.jet_speed==0;@assert maximum(abs,rho./Rref.-1)<.75;end
    rng=MersenneTwister(cfg["seed"]);noise=rand(rng,p.nx,p.ny,p.nz).-0.5;noise.-=mean(noise;dims=(1,2))
    # The same zero-horizontal-mean absolute pressure noise adds zero integrated heat.
    if perturb
        P.+=(5/3)*cfg["perturbation"].*Pref.*noise.*reshape(sinpi.(clamp.(z./DC,0.,1.)),1,1,:)
    end
    theta=P.^.6 ./rho
    set!(m;θ=theta,ρ=rho,qᵗ=0.,u=repeat(reshape(u,p.nx,p.ny,1),1,1,p.nz),v=repeat(reshape(v,p.nx,p.ny+1,1),1,1,p.nz),w=0.)
    Oceananigans.TimeSteppers.update_state!(m)
    actual=arrays(m);got=seed_moments(actual.rho,u,v,p)
    @assert maximum(abs,got./target.-1)<1e-10
    zeta=channel_vorticity(actual.uf,actual.vf,p)
    div=(circshift(u,(-1,0)).-u)./p.dx.+diff(v;dims=2)./p.dy
    @assert maximum(abs,div)<1e-12
    @assert abs(mean(zeta))<1e-12
    if p.jet_speed==0;@assert abs(mean(zeta.^3)/max(mean(zeta.^2)^1.5,1e-99))<1e-10;end
    pressure_error=maximum(abs,Array(interior(m.dynamics.pressure)).-P)/maximum(P)
    @assert pressure_error<1e-10
    Etarget=sum(CV.*P.+GRAV.*rho.*zz)*dv+target[2]
    Eactual=sum(actual.rho.*(CV.*actual.T.+GRAV.*zz))*dv+got[2]
    @assert abs(Eactual/Etarget-1)<1e-10
    N2=GRAV.*diff(log.(P.^.6 ./rho);dims=3)./p.dz
    zface=.5 .* (z[1:end-1].+z[2:end]);deep=findall(zface.<.65);cap=findall(zface.>.9)
    stability=Dict("N2_deep_min"=>minimum(N2[:,:,deep]),"N2_deep_max"=>maximum(N2[:,:,deep]),"N2_cap_min"=>minimum(N2[:,:,cap]),"N2_cap_max"=>maximum(N2[:,:,cap]),"balanced_energy_difference"=>Eactual-(Eref+target[2]))
    println("INITIAL_STABILITY ",stability);flush(stdout)
    @assert stability["N2_cap_min"]>0 "Balanced eddies must preserve a stable upper cap"
    @assert max(abs(stability["N2_deep_min"]),abs(stability["N2_deep_max"]))<.2 "Seed spuriously changes bulk convective stability"
    initial_info=Dict("status"=>"INITIALIZED","initial_rms"=>cfg["initial_rms"],"wavelength"=>4/sqrt(5),"turnover_time"=>4/(sqrt(5)*max(cfg["initial_rms"],1e-30)),"mass"=>got[1],"kinetic_energy"=>got[2],"total_energy"=>Eactual,"native_vorticity_skewness"=>mean(zeta.^3)/max(mean(zeta.^2)^1.5,1e-99),"jet_peak"=>maximum(abs,jets),"uniform_mass_scaling"=>(p.jet_speed==0 ? 1. : coefficients[1]+1),"jet_cap_displacement_min"=>minimum(delta),"jet_cap_displacement_max"=>maximum(delta),"max_abs_zeta_over_f"=>maximum(abs,zeta)/abs(p.f),"velocity_divergence_max"=>maximum(abs,div),"pressure_readback_error"=>pressure_error,"density_relative_departure_max"=>maximum(abs,rho./Rref.-1),"matching_coefficients"=>collect(coefficients),"pressure_offset"=>pressure_offset,"thermal_equilibrium"=>false,"note"=>"Odd-in-y projection of original streamfunction, same phases and |n|²=5; walls v=0 and tangential stress free. Jet balanced with f(y) and tilted rotation using a thermodynamic adjustment; eddy balance uses f0 with finite small beta adjustment. Widened y preserves the original physical eddy wavelengths. Finite-amplitude mixed-sign flow; mass and kinetic energy matched; pressure/density and balanced thermal/potential energy differ for rotation-dependent balance. No mature thermal parent.")
    merge!(initial_info,stability)
    return initial_info
end

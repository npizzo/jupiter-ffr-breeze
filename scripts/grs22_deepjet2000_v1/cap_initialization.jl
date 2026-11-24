# Solve the compressible balance characteristics for an x-uniform U(y,z).
# P_y = a P_z, a=fU/(g-hU); trace dz/dy=-a back to y=0.
function characteristic(y,z,p)
    n=max(1,ceil(Int,abs(y)/(p.dy/4)));ds=-y/n;v=z;D=1.;yy=y
    function rhs(y,z,D)
        J=jet_profile(y,p);U=jet_target(y,z,p);Uz=p.cap_jet ? J*cap_weight_z(z,p) : 0.
        f=p.f+p.beta*y;den=GRAV-p.fh*U
        (-f*U/den,-f*GRAV*Uz/den^2*D)
    end
    for _=1:n
        a,b=rhs(yy,v,D);c,d=rhs(yy+ds/2,v+ds*a/2,D+ds*b/2)
        e,f=rhs(yy+ds/2,v+ds*c/2,D+ds*d/2);g,h=rhs(yy+ds,v+ds*e,D+ds*f)
        v+=ds*(a+2*c+2*e+g)/6;D+=ds*(b+2*d+2*f+h)/6;yy+=ds
    end
    v,D
end
function apply_seed!(m,p,cfg;perturb=true)

    y=[(j-.5)*p.dy-p.ly/2 for j=1:p.ny];z=[zcell(k,p) for k=1:p.nz]
    q=[characteristic(yy,zz,p) for yy in y,zz in z]
    u=[jet_target(yy,zz,p) for yy in y,zz in z]
    P=[reference(q[j,k][1]).p for j=1:p.ny,k=1:p.nz]
    rho=[reference(q[j,k][1]).rho*GRAV*q[j,k][2]/(GRAV-p.fh*u[j,k]) for j=1:p.ny,k=1:p.nz]
    R=Array(interior(m.dynamics.reference_state.density));Pr=Array(interior(m.dynamics.reference_state.pressure))
    scale=mean(R)/mean(rho);P.*=scale;rho.*=scale
    # Enforce exactly the same native vertical pressure/reference-gravity balance.
    for k=p.nz-1:-1:1,j=1:p.ny
        rz=.5*(rho[j,k+1]+rho[j,k]);rref=.5*(R[1,1,k+1]+R[1,1,k])
        mz=.5*(rho[j,k+1]*u[j,k+1]+rho[j,k]*u[j,k])
        dp=Pr[1,1,k+1]-Pr[1,1,k]+p.dz*(-GRAV*(rz-rref)+p.fh*mz)
        P[j,k]=P[j,k+1]-dp
    end
    @assert minimum(rho)>0 && minimum(P)>0
    rep(a)=repeat(reshape(a,1,p.ny,p.nz),p.nx,1,1)
    set!(m;ρ=rep(rho),θ=rep(P.^.6 ./rho),u=rep(u),v=0.,w=0.,qᵗ=0.)
    Oceananigans.TimeSteppers.update_state!(m)
    if cfg["initial_rms"]>0;add_balanced_eddies!(m,p,cfg);end
    Dict{String,Any}("initial_rms"=>cfg["initial_rms"],"design"=>"Jet balanced by compressible characteristics and native vertical reintegration; full depth when jet_design=deep","uniform_mass_scaling"=>scale,"cap_low"=>p.cap_low,"cap_high"=>p.cap_high,"thermal_equilibrium"=>false)
end
function spatial_screen(m,p)
    a=arrays(m);r=dropdims(mean(a.rho;dims=1);dims=1);T=dropdims(mean(a.T;dims=1);dims=1);u=dropdims(mean(a.uf;dims=1);dims=1)
    P=dropdims(mean(Array(interior(m.dynamics.pressure));dims=1);dims=1)
    Pr=vec(Array(interior(m.dynamics.reference_state.pressure)));rr=vec(Array(interior(m.dynamics.reference_state.density)))
    l=log.((r.*T).^.6 ./r);ny,nz=size(r)
    uy=(u[3:end,2:end-1].-u[1:end-2,2:end-1])./(2p.dy)
    uz=(u[2:end-1,3:end].-u[2:end-1,1:end-2])./(2p.dz)
    sy=(l[3:end,2:end-1].-l[1:end-2,2:end-1])./(2p.dy)
    sz=(l[2:end-1,3:end].-l[2:end-1,1:end-2])./(2p.dz)
    f=reshape([p.f+p.beta*((j-.5)*p.dy-p.ly/2) for j=2:ny-1],ny-2,1)
    pv=(p.fh.+uz).*sy.+(f.-uy).*sz
    N2=GRAV.*sz;zz=[zcell(k,p) for k=2:nz-1];cap=findall(zz.>p.cap_low);deep=findall(zz.<.65)
    fy=reshape([p.f+p.beta*((j-1)*p.dy-p.ly/2) for j=2:ny],ny-1,1)
    ry=.5 .*(r[2:end,:].+r[1:end-1,:]);my=.5 .*((r.*u)[2:end,:].+(r.*u)[1:end-1,:])
    ay=.-(diff(P;dims=1)./p.dy.+fy.*my)./ry
    rz=.5 .*(r[:,2:end].+r[:,1:end-1]);mz=.5 .*((r.*u)[:,2:end].+(r.*u)[:,1:end-1])
    if p.neutral_coriolis
        fc=reshape([p.f+p.beta*((j-.5)*p.dy-p.ly/2) for j=1:ny],ny,1)
        my=.5 .*((r.*u.*fc)[2:end,:].+(r.*u.*fc)[1:end-1,:])
        ay=.-(diff(P;dims=1)./p.dy.+my)./ry
    end
    az=(-diff(P.-reshape(Pr,1,nz);dims=2)./p.dz.-GRAV.*(rz.-reshape(.5 .*(rr[2:end].+rr[1:end-1]),1,nz-1)).+p.fh.*mz)./rz
    scale=abs(p.f)*sqrt(mean(u.^2)); normalized=scale>1e-10
    Dict("cap_PV_min"=>minimum(pv[:,cap]),"cap_PV_max"=>maximum(pv[:,cap]),"cap_N2_min"=>minimum(N2[:,cap]),"deep_PV_min"=>minimum(pv[:,deep]),"deep_PV_max"=>maximum(pv[:,deep]),"vertical_absolute_vorticity_min"=>minimum(f.-uy),"balance_normalization_defined"=>normalized,"horizontal_balance_rms"=>sqrt(mean(ay.^2)),"vertical_balance_max"=>maximum(abs,az),"horizontal_balance_rms_over_fU"=>(normalized ? sqrt(mean(ay.^2))/scale : 0.),"vertical_balance_max_over_fU"=>(normalized ? maximum(abs,az)/scale : 0.),"top_jet_dudz_max"=>maximum(abs,(u[:,end].-u[:,end-1])./p.dz),"bottom_jet_dudz_max"=>maximum(abs,(u[:,2].-u[:,1])./p.dz),"Tmin"=>minimum(T),"Tmax"=>maximum(T),"theta_min"=>minimum(exp.(l)),"theta_max"=>maximum(exp.(l)))
end

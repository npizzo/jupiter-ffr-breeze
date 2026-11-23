using Oceananigans.BoundaryConditions: fill_halo_regions!
function source_checks(cfg)
    c=merge(cfg,Dict("nx"=>8,"ny"=>16,"nz"=>128,"lx"=>.25,"initial_rms"=>0.,"jet_speed"=>0.,"perturbation"=>0.))
    m,p=build_model(c;arch=CPU(),perturb=false)
    T=[reference(zcell(k,p)).T+.03cospi((j-.5)/p.ny) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    r=[reference(zcell(k,p)).rho for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    set!(p.nu,(x,y,z)->.001*(1+.2cospi(y/4)));fill_halo_regions!(p.nu)
    F=(T=T,ρᵈ=r);errors=Dict{String,Float64}()
    # Each external source independently telescopes to its stated boundary/volume power.
    for (name,b,c,t,cap) in [("off",0.,0.,0.,0.),("bottom",1.,0.,0.,0.),("interior_conduction",0.,1.,0.,0.),("top",0.,0.,1.,0.),("cap",0.,0.,0.,1.),("all",1.,1.,1.,1.)]
        pp=merge(p,(;bottom_factor=b,conduction_factor=c,top_factor=t,cap_factor=cap))
        Q=[external_heating(i,j,k,F,pp) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
        expected=p.lx*p.ly*(b*p.flux-t*mean(2p.kc.*(T[:,:,end].-p.tt)./p.dz))-cap*sum(r.*(T.-p.tt).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,p.nz))*p.dx*p.dy*p.dz
        errors[name]=abs(sum(Q)*p.dx*p.dy*p.dz-expected);@assert errors[name]<1e-10
    end
    Q=[sgs_heating(i,j,k,F,p) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
    errors["closed_entropy_flux"]=abs(sum(Q)*p.dx*p.dy*p.dz);@assert errors["closed_entropy_flux"]<1e-10
    @assert parameters(merge(c,Dict("bottom_flux"=>0.))).kc==p.kc
    # Volumetric heating conversion: d(cv*P)/d(rho*theta) = cp*T/theta.
    theta=(r.*T).^.6 ./r;testQ=ones(size(T));tendency=theta.*testQ./(2.5 .*T)
    errors["heating_chain_rule"]=maximum(abs,2.5 .*T./theta.*tendency.-testQ);@assert errors["heating_chain_rule"]<1e-14
    Dict("status"=>"PASS","errors"=>errors,"conductivity_independent_of_bottom_flux"=>true)
end

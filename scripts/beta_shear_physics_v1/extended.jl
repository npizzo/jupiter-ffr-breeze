using Oceananigans.Solvers: FFTW
using .FFTW: fft, dct
function extended_output(m,p,cfg,out,t)
    a=arrays(m);r=a.rho;T=a.T;nu=Array(interior(p.nu))
    rt=Array(interior(m.formulation.potential_temperature_density))
    dye=Array(interior(m.tracers.dye))./r
    path=joinpath(out,"extrema.csv")
    if !isfile(path)
        open(path,"w") do io;println(io,"time,rhotheta_min,rhotheta_max,theta_min,theta_max,tracer_min,tracer_max,nu_max,cap_nu_mean");end
    end
    cap=findall([zcell(k,p)>.82 for k=1:p.nz])
    open(path,"a") do io;println(io,join((t,minimum(rt),maximum(rt),minimum(rt./r),maximum(rt./r),minimum(dye),maximum(dye),maximum(nu),mean(nu[:,:,cap])),','));end
    if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
        # Horizontal fluxes telescope. This is the exact horizontal-mean SGS heat tendency.
        s=CV.*log.(T).-log.(r);c=r.*T.*nu./p.pr
        qz=zeros(p.nx,p.ny,p.nz+1)
        qz[:,:,2:end-1].=-.5 .*(c[:,:,1:end-1].+c[:,:,2:end]).*diff(s;dims=3)./p.dz
        entropy_profile=vec(mean(-diff(qz;dims=3)./p.dz;dims=(1,2))).*p.entropy_diffusion_factor
        n2=GRAV.*diff(log.(rt./r);dims=3)./p.dz
        # Directional spectra of sqrt(rho)*cell velocity. x is periodic; y uses orthonormal DCT-II.
        ex=zeros(p.nx);ey=zeros(p.ny);energy=0.
        for vel in (a.u,a.v,a.w)
            q=sqrt.(r).*vel;energy+=.5mean(q.^2)
            fx=fft(q,1);ex.+=.5vec(sum(abs2.(fx);dims=(2,3)))./(p.nx^2*p.ny*p.nz)
            fy=dct(q,2);ey.+=.5vec(sum(abs2.(fy);dims=(1,3)))./(p.nx*p.ny*p.nz)
        end
        @assert abs(sum(ex)-energy)<1e-10*max(energy,1.)
        @assert abs(sum(ey)-energy)<1e-10*max(energy,1.)
        jldsave(joinpath(out,@sprintf("profiles_%07d.jld2",round(Int,100t)));t,z=[zcell(k,p) for k=1:p.nz],
            mean_temperature=vec(mean(T;dims=(1,2))),mean_density=vec(mean(r;dims=(1,2))),
            mean_u=dropdims(mean(a.uf;dims=1);dims=1),mean_nu=vec(mean(nu;dims=(1,2))),
            mean_N2=vec(mean(n2;dims=(1,2))),entropy_profile,spectrum_x=ex,spectrum_y=ey,
            spectrum_energy=energy,nu_max=maximum(nu))
    end
end

function implementation_checks(base)
    c=merge(base,Dict("nx"=>8,"ny"=>32,"nz"=>80,"lx"=>.25,"initial_rms"=>0.,"jet_speed"=>.015625,"jet_design"=>"cap","neutral_coriolis"=>true))
    m,p=build_model(c;arch=CPU(),perturb=false)
    # A strongly variable density and all velocity components exercise the weighted pairing.
    set!(m;u=(x,y,z)->.03cospi(2x/p.lx)*cospi(y/p.ly),v=(x,y,z)->.02sinpi(y/p.ly)*cospi(2x/p.lx),w=(x,y,z)->.01sinpi(z/p.lz)*cospi(y/p.ly))
    Oceananigans.TimeSteppers.update_state!(m)
    a=audit(m,p);@assert abs(a["coriolis_work"])<1e-11
    # Verify the native closure receives g*d(log theta)/dz, with consistent stencil.
    B=Breeze.AtmosphereModels
    b=Oceananigans.TurbulenceClosures.buoyancy_force(m)
    tr=Oceananigans.TurbulenceClosures.buoyancy_tracers(m)
    f=Oceananigans.BuoyancyFormulations.∂z_b
    theta=Array(interior(m.formulation.potential_temperature))
    err=maximum(abs(f(2,j,k,m.grid,b,tr)-GRAV*(log(theta[2,j,k])-log(theta[2,j,k-1]))/p.dz) for j=2:p.ny-1,k=2:p.nz)
    @assert err<1e-8
    # Verify native Lilly multiplier, including the installed Sigma^2 convention.
    l,lp=build_model(merge(c,Dict("closure"=>"lilly"));arch=CPU(),perturb=false)
    o,op=build_model(merge(c,Dict("closure"=>"original"));arch=CPU(),perturb=false)
    nuc=Array(interior(l.closure_fields.νₑ));nuo=Array(interior(o.closure_fields.νₑ))
    cap=findall([zcell(k,lp)>.85 for k=1:lp.nz])
    @assert maximum(nuo[:,:,cap])>0
    @assert maximum(nuc[:,:,cap])<1e-12
    Dict("status"=>"PASS","coriolis_power"=>a["coriolis_work"],"native_N2_error"=>err,
         "original_cap_nu_max"=>maximum(nuo[:,:,cap]),"lilly_cap_nu_max"=>maximum(nuc[:,:,cap]),
         "lilly_definition"=>"nu=(Cs*Delta)^2*sqrt(2*S:S)*sqrt(max(0,1-N2_positive/(Pr*2*S:S))); installed Cb=0.5/Pr")
end

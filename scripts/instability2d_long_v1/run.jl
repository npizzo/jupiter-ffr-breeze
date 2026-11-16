include("audit.jl")
include("spectra.jl")
using JLD2, SHA
out,case,source=ARGS
case in ("forced","unforced")||error("Unknown case")
ispath(out)&&error("Refusing to overwrite existing output");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];t0=d["t"]
n=128;dt=0.01;sigma=0.6;target=1200.;steps=round(Int,(target/sigma-t0)/dt)
o=BandOperator(n,n;kmin=6.,kmax=8.,allow_zonal=false)
m=build_model(o;f=1.);seed!(m,u0,v0);m.clock.time=t0
b=Audit(o;a=case=="forced" ? 0.009375 : 0.,dt)
config=Dict("case"=>case,"grid"=>[n,n],"domain"=>[2pi,2pi],"f"=>1.,"planetary_beta"=>0.,
    "dt"=>dt,"target_growth_times"=>target,"branch_growth_time"=>sigma*t0,
    "target_iterations"=>steps,"amplitude"=>b.amplitude,"viscosity"=>b.nu,"cubic_velocity_drag"=>b.alpha,
    "forcing_band"=>[6.,8.],"allow_zonal_forcing"=>false,"remove_mean_drag"=>true,
    "initial_rms_velocity"=>sqrt(mean(u0.^2 .+v0.^2)),"sigma_reference"=>sigma,
    "source_checkpoint"=>source,"source_checkpoint_sha256"=>bytes2hex(sha256(read(source))),
    "physical_time_mapping"=>"none; sigma*t uses the same reference sigma=0.6 in BOTH cases, including unforced",
    "compute_location"=>"Unity","wall_time_limit_seconds"=>12600,
    "spectral_flux_sign"=>"Pi(K)=-sum(k<=K,T_advection); negative means net upscale transfer across K",
    "spectral_shells"=>"integer upper edges: K-1 < |k| <= K",
    "spectral_cadence_iterations"=>100,"frame_cadence_iterations"=>200,
    "oceananigans_version"=>string(pkgversion(Oceananigans)),
    "source_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in [@__FILE__,joinpath(@__DIR__,"audit.jl"),joinpath(@__DIR__,"spectra.jl"),joinpath(@__DIR__,"../instability2d_v1/core.jl")]))
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"statistics.csv"),"w") do io
    println(io,"iteration,t,growth_times,energy,enstrophy,skewness,zeta_min,zeta_max,low_k_energy_fraction,integral_length,energy_advection,energy_injection,energy_viscosity,energy_drag,numerical_advection_M3_fraction")
end
open(joinpath(out,"spectra.csv"),"w") do io
    println(io,"iteration,t,growth_times,k,energy,T_advection,T_injection,T_viscosity,T_drag,Pi_advection")
end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,growth_times,file");end
started=time();lastframe=Ref(-1);laststat=Ref(-1);maxparseval=Ref(0.);maxwork=Ref(0.)
function audited_forcing!(model)
    b(model)
    if model.clock.stage==1 && model.clock.iteration%100==0
        st=b.stages[1];u,v=state(model)
        s=spectral_checks(o,st.q,st.parts,u,v)
        maxparseval[]=max(maxparseval[],s.eerr);maxwork[]=max(maxwork[],s.terr)
        write_spectra(joinpath(out,"spectra.csv"),model.clock,s.s)
    end
end
function output!(sim;force=false)
    m=sim.model;finish!(b,m);it=m.clock.iteration
    if (it%100==0 || force) && laststat[]!=it
        u,v=state(m);q=vorticity(o,u,v);qr=sqrt(mean(q.^2));energy=mean(u.^2 .+v.^2)/2
        qh=FFTW.fft(q);invlam=map(x->x==0 ? 0.0 : 1/x,o.lambda)
        ek=abs2.(qh).*invlam;low=sum(ek[o.k2.<=4])/max(sum(ek),1e-30)
        invk=map(x->x==0 ? 0.0 : 1/sqrt(x),o.k2);lengthscale=2pi*sum(ek.*invk)/max(sum(ek),1e-30)
        ratio=b.absolute_third[1]/max(sum(b.absolute_third[2:4]),1e-30)
        open(joinpath(out,"statistics.csv"),"a") do io
            println(io,join((it,m.clock.time,sigma*m.clock.time,energy,qr^2/2,mean(q.^3)/max(qr^3,1e-30),minimum(q),maximum(q),low,lengthscale,b.energy...,ratio),','))
        end
        laststat[]=it
    end
    if (it%200==0 || force) && lastframe[]!=it
        u,v=state(m);q=vorticity(o,u,v);name=@sprintf("zeta_%07d.f32",it)
        open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
        open(joinpath(out,"frames.csv"),"a") do io;println(io,join((it,m.clock.time,sigma*m.clock.time,name),','));end
        lastframe[]=it
    end
    if it>0 && it%2000==0
        u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),iteration=it,config)
        mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true)
        println(@sprintf("PROGRESS case=%s growth_times=%.2f elapsed=%.1fs CFL=%.3f energy_closure=%.3g",case,sigma*m.clock.time,time()-started,b.max_cfl,b.max_energy_error));flush(stdout)
    end
end
sim=Simulation(m;Δt=dt,stop_iteration=steps,wall_time_limit=12600.,verbose=false)
add_callback!(sim,audited_forcing!,IterationInterval(1);name=:audited_forcing,callsite=TendencyCallsite())
add_callback!(sim,output!,IterationInterval(1);name=:paired_output)
output!(sim;force=true);run!(sim);output!(sim;force=true)
u,v=state(m);jldsave(joinpath(out,"checkpoint.jld2");u,v,t=Float64(m.clock.time),iteration=m.clock.iteration,config)
result=report(b);result["status"]=m.clock.iteration>=steps ? "COMPLETED" : "PARTIAL_WALL_TIME_LIMIT"
result["growth_times_completed"]=sigma*m.clock.time;result["elapsed_seconds"]=time()-started
result["max_spectral_energy_error"]=maxparseval[];result["max_spectral_work_error"]=maxwork[]
result["scope"]="Paired mature-state forcing switch; broad inertial-range or sign-selection claims require additional evidence"
open(joinpath(out,"review.toml"),"w") do io;TOML.print(io,result;sorted=true);end
println("LONG_PAIR_FINISHED");TOML.print(stdout,result;sorted=true)
result["status"]=="COMPLETED"||error("Partial output retained; no automatic extension")

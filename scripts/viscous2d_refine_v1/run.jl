include("../stirring2d_v1/forcing.jl")
using JLD2, SHA
function write_native_spectra(path,clock,s)
    open(path,"a") do io
        for j in eachindex(s.k)
            println(io,join((clock.iteration,clock.time,clock.time,s.k[j],s.E[j],s.T[j,:]...,s.flux[j]),','))
        end
    end
end
out,case,source=ARGS[1:3]
probe=length(ARGS)==4 && ARGS[4]=="probe"
scale=0.0002/0.009375
case=="unforced"||error("Only the authorized unforced check is supported")
ispath(out)&&error("Refusing to overwrite existing output");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"];t0=d["t"]
n=512;dt=0.0025;sigma=1.0;target=8.;steps=round(Int,(target/sigma-t0)/dt)
o=BandOperator(n,n;kmin=12.,kmax=16.,allow_zonal=false)
m=build_model(o;f=1.);seed!(m,u0,v0);m.clock.time=t0
b=Audit(o;a=case=="forced" ? 0.0002 : 0.,nu=0.0001,alpha=0.,dt)
epsilon=scale*reference_injection(o,u0,v0)
stir=Stirring(o;case=case in ("forced","unforced") ? "random_sign" : case,epsilon)
stochastic=case in ("random_sign","cyclonic_sign")
kicksteps=round(Int,stir.interval/dt);@assert isapprox(kicksteps*dt,stir.interval)
E0=mean(u0.^2 .+v0.^2)/2;M30=mean(vorticity(o,u0,v0).^3)
config=Dict("case"=>case,"grid"=>[n,n],"domain"=>[2pi,2pi],"f"=>1.,"planetary_beta"=>0.,
    "dt"=>dt,"target_model_time"=>target,"initial_model_time"=>sigma*t0,
    "target_iterations"=>steps,"amplitude"=>b.amplitude,"viscosity"=>b.nu,"cubic_velocity_drag"=>b.alpha,
    "forcing_band"=>[12.,16.],"allow_zonal_forcing"=>false,"remove_mean_drag"=>true,
    "initial_rms_velocity"=>sqrt(mean(u0.^2 .+v0.^2)),"sigma_reference"=>sigma,
    "source_checkpoint"=>source,"source_checkpoint_sha256"=>bytes2hex(sha256(read(source))),
    "physical_time_mapping"=>"Dimensionless model time; no mapping to planetary days",
    "compute_location"=>"Unity","wall_time_limit_seconds"=>5400,
    "spectral_flux_sign"=>"Pi(K)=-sum(k<=K,T_advection); negative means net upscale transfer across K",
    "spectral_shells"=>"integer upper edges: K-1 < |k| <= K",
    "spectral_cadence_iterations"=>81,"frame_cadence_iterations"=>200,
    "oceananigans_version"=>string(pkgversion(Oceananigans)),
    "stirring"=>(stochastic ? "independent fixed-amplitude band-limited vortex pulses; compensated circulation" : case=="forced" ? "continuous instability forcing" : "unforced decay"),
    "forcing_scale"=>scale,"diagnostic_phase_cycle_samples"=>80,
    "diagnostic_sampling"=>"Every 81 steps at dt=.0025; all 80 discrete pulse phases sampled equally over complete 80-sample cycles",
    "movie_sampling"=>"Unforced snapshots every 0.5 model time",
    "probe"=>probe,
    "stochastic_pulses_enabled"=>stochastic,"pulse_interval"=>stir.interval,"pulse_amplitude"=>(stochastic ? stir.amplitude : 0.),"expected_energy_injection_rate"=>(case=="unforced" ? 0. : epsilon),
    "energy_match"=>"weak a=.0002 instability work at identical initial random state; expectation matched, realized cross work measured separately",
    "position_seed"=>23017,"sign_seed"=>23018,"sign_rule"=>(case in ("forced","unforced") ? "no imposed polarity" : case=="random_sign" ? "independent equiprobable +/- cores" : "core sign = sign(f)"),
    "pulse_template_skewness"=>mean(stir.q.^3)/mean(stir.q.^2)^1.5,
    "expected_third_moment_injection_rate"=>(case!="cyclonic_sign" ? 0.0 : sign(stir.f)*stir.amplitude^3*mean(stir.q.^3)/stir.interval),
    "source_sha256"=>Dict(relpath(p,@__DIR__)=>bytes2hex(sha256(read(p))) for p in [@__FILE__,joinpath(@__DIR__,"../stirring2d_v1/forcing.jl"),joinpath(@__DIR__,"../instability2d_long_v1/audit.jl"),joinpath(@__DIR__,"../instability2d_long_v1/spectra.jl"),joinpath(@__DIR__,"../instability2d_v1/core.jl")]))
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"statistics.csv"),"w") do io
    println(io,"iteration,t,model_time,energy,enstrophy,skewness,zeta_min,zeta_max,low_k_energy_fraction,integral_length,energy_advection,energy_injection,energy_viscosity,energy_drag,numerical_advection_M3_fraction,pulse_phase_step,lowpass_skewness,urms,eddy_turnover_time,reference_energy_input_time,pulse_velocity_ratio,rms_native_vorticity_tendency,rms_continuous_forcing_tendency")
end
open(joinpath(out,"spectra.csv"),"w") do io
    println(io,"iteration,t,model_time,k,energy,T_advection,T_injection,T_viscosity,T_drag,Pi_advection")
end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,model_time,file");end
open(joinpath(out,"events.csv"),"w") do io
    println(io,"event,iteration,t,model_time,shift_x,shift_y,core_sign,energy_quadratic,energy_cross,energy_total,third_moment_increment")
end
open(joinpath(out,"kick_spectra.csv"),"w") do io;println(io,"event,iteration,t,model_time,k,energy_increment");end
started=time();lastframe=Ref(-1);laststat=Ref(-1);lastspec=Ref(-1);maxparseval=Ref(0.);maxwork=Ref(0.)
function audited_forcing!(model)
    b(model)
end
function output!(sim;force=false)
    m=sim.model;finish!(b,m);it=m.clock.iteration
    if stochastic && it>0 && it%kicksteps==0 && stir.last_iteration!=it
        ev=kick!(stir,b,m)
        open(joinpath(out,"events.csv"),"a") do io
            println(io,join((stir.events,it,m.clock.time,sigma*m.clock.time,ev.ij...,ev.sg,ev.quadratic,ev.cross,ev.de,ev.dm3),','))
        end
        open(joinpath(out,"kick_spectra.csv"),"a") do io
            for k=0:length(ev.shellenergy)-1;println(io,join((stir.events,it,m.clock.time,sigma*m.clock.time,k,ev.shellenergy[k+1]),','));end
        end
    end
    if it%81==0 && lastspec[]!=it && b.iteration==it
        st=b.stages[1];u,v=state(m);s=spectral_checks(o,st.q,st.parts,u,v)
        maxparseval[]=max(maxparseval[],s.eerr);maxwork[]=max(maxwork[],s.terr)
        write_native_spectra(joinpath(out,"spectra.csv"),m.clock,s.s);lastspec[]=it
    end
    if (it%81==0 || force) && laststat[]!=it
        u,v=state(m);q=vorticity(o,u,v);qr=sqrt(mean(q.^2));energy=mean(u.^2 .+v.^2)/2
        qh=FFTW.fft(q);invlam=map(x->x==0 ? 0.0 : 1/x,o.lambda)
        ek=abs2.(qh).*invlam;low=sum(ek[o.k2.<=4])/max(sum(ek),1e-30)
        invk=map(x->x==0 ? 0.0 : 1/sqrt(x),o.k2);lengthscale=2pi*sum(ek.*invk)/max(sum(ek),1e-30)
        ql=spectral_curl(o,q,(o.k2.>0).&(o.k2.<=4));sl=mean(ql.^3)/max(mean(ql.^2)^1.5,1e-30)
        urms=sqrt(2energy);tau=lengthscale/max(urms,1e-30)
        native_rms=b.iteration==it ? sqrt(mean(b.stages[1].parts[1].^2)) : 0.
        forcing_rms=sqrt(mean(band_curl(o,q,b.amplitude).^2))
        ratio=b.absolute_third[1]/max(sum(b.absolute_third[2:4]),1e-30)
        abs(energy-E0-sum(b.energy))/E0<1e-8||error("Total RK plus pulse energy did not close")
        abs(mean(q.^3)-M30-sum(b.third))/max(mean(abs.(q).^3),1e-10)<1e-7||error("Total RK plus pulse third moment did not close")
        open(joinpath(out,"statistics.csv"),"a") do io
            println(io,join((it,m.clock.time,sigma*m.clock.time,energy,qr^2/2,mean(q.^3)/max(qr^3,1e-30),minimum(q),maximum(q),low,lengthscale,b.energy...,ratio,it%kicksteps,sl,urms,tau,energy/epsilon,stochastic ? stir.amplitude/max(urms,1e-30) : 0.,native_rms,forcing_rms),','))
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
        u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),iteration=it,config,positions_rng=stir.positions,signs_rng=stir.signs,event_count=stir.events)
        mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true)
        println(@sprintf("PROGRESS case=%s model_time=%.2f elapsed=%.1fs CFL=%.3f energy_closure=%.3g",case,sigma*m.clock.time,time()-started,b.max_cfl,b.max_energy_error));flush(stdout)
    end
end
sim=Simulation(m;Δt=dt,stop_iteration=steps,wall_time_limit=5400.,verbose=false)
add_callback!(sim,audited_forcing!,IterationInterval(1);name=:audited_forcing,callsite=TendencyCallsite())
add_callback!(sim,output!,IterationInterval(1);name=:paired_output)
output!(sim;force=true);run!(sim);output!(sim;force=true)
u,v=state(m);jldsave(joinpath(out,"checkpoint.jld2");u,v,t=Float64(m.clock.time),iteration=m.clock.iteration,config,positions_rng=stir.positions,signs_rng=stir.signs,event_count=stir.events)
result=report(b);result["status"]=m.clock.iteration>=steps ? "COMPLETED" : "PARTIAL_WALL_TIME_LIMIT"
result["energy_terms"]=["advection_coriolis",stochastic ? "stochastic_pulses" : "instability","viscosity","cubic_drag"]
result["third_moment_terms"]=result["energy_terms"]
result["min_interpulse_injection_work"]=pop!(result,"min_injection_work")
result["model_time_completed"]=sigma*m.clock.time;result["elapsed_seconds"]=time()-started
result["max_spectral_energy_error"]=maxparseval[];result["max_spectral_work_error"]=maxwork[]
result["pulse_events"]=stir.events;result["positive_core_events"]=stir.positive_events
result["pulse_quadratic_energy"]=stir.energy_quadratic;result["pulse_cross_energy"]=stir.energy_cross
result["realized_mean_energy_injection"]=b.energy[2]/(m.clock.time-t0)
result["reference_mean_energy_injection"]=epsilon
result["forcing_scale"]=scale
result["scope"]="512-square unforced spatial refinement; same dt, viscosity and initial vorticity Fourier modes as the 256-square check"
open(joinpath(out,"review.toml"),"w") do io;TOML.print(io,result;sorted=true);end
println("LONG_PAIR_FINISHED");TOML.print(stdout,result;sorted=true)
result["status"]=="COMPLETED"||error("Partial output retained; no automatic extension")

include("forcing.jl")
using Test,JLD2,SHA
out,source=ARGS;ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"]
o=BandOperator(128,128;kmin=6.,kmax=8.,allow_zonal=false)
epsilon=reference_injection(o,u0,v0);metrics=Dict{String,Any}()
function evolve(case;dt=0.005,steps=4000)
    m=build_model(o;f=1.);seed!(m,u0,v0);b=Audit(o;a=0.,dt)
    s=Stirring(o;case,epsilon);ks=round(Int,s.interval/dt)
    E0=mean(u0.^2 .+v0.^2)/2;M0=mean(vorticity(o,u0,v0).^3)
    maxerr=Ref(0.);maxmerr=Ref(0.);locations=Tuple{Int,Int}[];signs=Int[]
    function finish_and_kick!(sim)
        finish!(b,sim.model);it=sim.model.clock.iteration
        if it>0 && it%ks==0 && s.last_iteration!=it
            ev=kick!(s,b,sim.model);push!(locations,ev.ij);push!(signs,ev.sg)
        end
        if it>0 && it%100==0
            u,v=state(sim.model);q=vorticity(o,u,v);en=mean(u.^2 .+v.^2)/2
            maxerr[]=max(maxerr[],abs(en-E0-sum(b.energy))/E0)
            maxmerr[]=max(maxmerr[],abs(mean(q.^3)-M0-sum(b.third))/max(mean(abs.(q).^3),1e-10))
        end
    end
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:audit,callsite=TendencyCallsite())
    add_callback!(sim,finish_and_kick!,IterationInterval(1);name=:kicks)
    run!(sim);finish_and_kick!(sim)
    (;m,b,s,locations,signs,energy_error=maxerr[],third_error=maxmerr[])
end
@testset "Stochastic polarity controls" begin
    s=Stirring(o;case="random_sign",epsilon)
    @test epsilon>0
    @test abs(mean(s.q))<1e-12
    @test abs(mean(s.u))+abs(mean(s.v))<1e-12
    @test norm(vorticity(o,s.u,s.v)-s.q)/norm(s.q)<1e-12
    div=(circshift(s.u,(-1,0)).-s.u)./(o.lx/o.nx).+(circshift(s.v,(0,-1)).-s.v)./(o.ly/o.ny)
    @test maximum(abs,div)<1e-12
    @test s.q[1,1]==maximum(s.q)
    @test minimum(s.q)<0
    @test mean(s.q.^3)>0
    @test abs(mean(s.u.^2 .+s.v.^2)-1)<1e-12
    @test abs(s.amplitude^2/2-epsilon*s.interval)<1e-12
    @test maximum(abs,FFTW.fft(s.q)[.!o.mask])<1e-9
    # All possible grid translations have exactly zero mean cross-work,
    # proving equal EXPECTED injection for either sign law, independently of u.
    crossmap=real.(FFTW.ifft(conj.(FFTW.fft(u0)).*FFTW.fft(s.amplitude.*s.u).+
                           conj.(FFTW.fft(v0)).*FFTW.fft(s.amplitude.*s.v)))./length(u0)
    @test abs(mean(crossmap))<1e-14
    # Pointwise f/sign covariance and full velocity reflection at a chosen center.
    sm=Stirring(o;case="cyclonic_sign",epsilon,f=-1.)
    ij=(13,19);du=circshift(s.u,ij);dv=circshift(s.v,ij)
    ur,vr=reflect_state(du,dv)
    @test norm(ur+circshift(sm.u,(-ij[1],ij[2])))/norm(ur)<1e-12
    @test norm(vr+circshift(sm.v,(-ij[1],ij[2])))/norm(vr)<1e-12
    # Same random locations, independent fair coin; fixed case always f-aligned.
    sr=Stirring(o;case="random_sign",epsilon);sp=Stirring(o;case="cyclonic_sign",epsilon)
    signs=Int[]
    for i=1:256
        ijr,sg=draw_event!(sr);ijp,sgp=draw_event!(sp);ijm,sgm=draw_event!(sm)
        @test ijr==ijp==ijm
        @test sgp==1 && sgm==-1
        push!(signs,sg)
    end
    @test -1 in signs && 1 in signs
    # Integrated 20-time-unit probes from the actual shared turbulent checkpoint.
    runs=Dict{String,Any}()
    for case in ("random_sign","cyclonic_sign")
        r=evolve(case);runs[case]=r
        @test r.s.events==100
        @test r.b.finished==3999
        @test r.energy_error<1e-8
        @test r.third_error<1e-7
        @test r.b.max_energy_error<1e-9
        @test r.b.max_increment_error<1e-9
        @test r.b.max_cfl<0.5
        metrics[case]=Dict("events"=>r.s.events,"positive_events"=>r.s.positive_events,
            "max_total_energy_error"=>r.energy_error,"max_total_third_moment_error"=>r.third_error,
            "max_cfl"=>r.b.max_cfl,"quadratic_energy"=>r.s.energy_quadratic,"cross_energy"=>r.s.energy_cross,
            "native_rk_energy_error"=>r.b.max_energy_error)
    end
    @test runs["random_sign"].locations==runs["cyclonic_sign"].locations
    @test runs["random_sign"].s.energy_quadratic==runs["cyclonic_sign"].s.energy_quadratic
    # Pulse times stay physical when halving dt; test the inter-pulse integrator.
    a=evolve("cyclonic_sign";steps=200);h=evolve("cyclonic_sign";dt=.0025,steps=400)
    u,v=state(a.m);uh,vh=state(h.m)
    err=sqrt(mean((u.-uh).^2 .+(v.-vh).^2))/sqrt(mean(uh.^2 .+vh.^2))
    @test a.locations==h.locations && a.s.events==h.s.events==5
    @test err<1e-3
    metrics["one_time_unit_dt_refinement_error"]=err
    metrics["expected_energy_injection_rate"]=epsilon
    metrics["pulse_template_skewness"]=mean(s.q.^3)/mean(s.q.^2)^1.5
    metrics["translation_mean_cross_work"]=mean(crossmap)
end
metrics["status"]="PASS";metrics["compute_location"]="Unity"
metrics["source_checkpoint_sha256"]=bytes2hex(sha256(read(source)))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
println("STIRRING_CHECK_PASS");TOML.print(stdout,metrics;sorted=true)

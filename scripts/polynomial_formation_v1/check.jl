include("model.jl")
using Test
initial,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(initial);u0,v0=d["u"],d["v"]
@assert d["t"]==0. && size(u0)==size(v0)==(512,512)
metrics=Dict{String,Any}();rel(a,b)=sqrt(sum(abs2,a.-b)/max(sum(abs2,b),eps()))
function evolve(u,v;steps=20,dt=.0025,t=0.)
    m,b=setup_noadvection(u,v;t,dt)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:polynomial,callsite=TendencyCallsite())
    run!(sim);(;m,b,sim)
end
@testset "Only-advection intervention and native reconstruction" begin
    errors=Float64[]
    # Existing states spanning formation: validate the Python flux algebra against native operators.
    for (case,t) in [("biased_plus",0.),("biased_plus",40.),("biased_plus",80.),("biased_plus",160.),("biased_plus",240.),("biased_plus",400.),("cubic_control",160.),("biased_minus",160.)]
        ending=max(80,80ceil(Int,t/80));p=joinpath("data/viscous2d-feedback512-v1-1925368",case,@sprintf("to_%03d",ending),@sprintf("zeta_%08d.f32",round(Int,2t)))
        q=Float64.(reshape(reinterpret(Float32,read(p)),512,512))
        o=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false);u,v=velocity_from_curl(o,q)
        m,b=setup_maintenance(u,v;case="full",t)
        dx=o.lx/o.nx;uc=.5 .*(u.+circshift(u,(-1,0)));vc=.5 .*(v.+circshift(v,(0,-1)))
        cross=.25 .*(u.+circshift(u,(0,1))).*(v.+circshift(v,(1,0)))
        au=-(uc.^2 .-circshift(uc.^2,(1,0)) .+circshift(cross,(0,-1)).-cross)./dx
        av=-(circshift(cross,(-1,0)).-cross .+vc.^2 .-circshift(vc.^2,(0,1)))./dx
        err=rel(native_advection_curl(m,o),vorticity(o,au,av));push!(errors,err);@test err<1e-9
    end
    m,b=setup_noadvection(u0,v0);o=b.base.op
    @test maximum(abs,native_advection_curl(m,o))==0.
    @test abs(mean(u0))+abs(mean(v0))<1e-12
    captured=Ref(false);native_error=Ref(NaN)
    sim=Simulation(m;Δt=.0025,stop_iteration=1,verbose=false)
    function capture(m)
        if !captured[]
            u,v=state(m);q=vorticity(o,u,v)
            gu=copy(dropdims(Array(interior(m.timestepper.Gⁿ.u));dims=3));gv=copy(dropdims(Array(interior(m.timestepper.Gⁿ.v));dims=3))
            native_error[]=rel(vorticity(o,gu,gv),spectral_curl(o,q,-.0001 .*o.lambda));captured[]=true
        end
        b(m)
    end
    add_callback!(sim,capture,IterationInterval(1);name=:capture,callsite=TendencyCallsite());run!(sim)
    @test captured[] && native_error[]<1e-8
    coarse=evolve(u0,v0);u,v=state(coarse.m)
    fine=evolve(u0,v0;steps=40,dt=.00125);uf,vf=state(fine.m)
    dt_error=max(rel(u,uf),rel(v,vf));@test dt_error<1e-5
    half=evolve(u0,v0;steps=9);uh,vh=state(half.m)
    restart_path=joinpath(out,"restart.jld2");jldsave(restart_path;u=uh,v=vh,t=Float64(half.m.clock.time))
    rd=load(restart_path);res=evolve(rd["u"],rd["v"];steps=11,t=rd["t"]);ur,vr=state(res.m)
    restart_error=max(rel(u,ur),rel(v,vr));@test restart_error<1e-11
    q0=vorticity(o,u0,v0);q1=vorticity(o,u,v)
    p0=potential_diagnostics(o,q0);p1=potential_diagnostics(o,q1)
    @test p1.potential<p0.potential
    src=maintenance_sources(o,q0);g=-(src.linear.+src.viscous.+src.quadratic.+src.cubic)
    @test isapprox(mean(g.^2),p0.gradient_norm_squared;rtol=1e-10)
    coarse.sim.stop_iteration=70;GC.gc();perf=@timed run!(coarse.sim)
    @test perf.time/50<.065 && perf.bytes/50<65536
    merge!(metrics,Dict("native_reconstruction_max_error"=>maximum(errors),"native_diffusion_only_error"=>native_error[],
        "restart_relative_error"=>restart_error,"half_dt_relative_error"=>dt_error,
        "seconds_per_step"=>perf.time/50,"bytes_per_step"=>perf.bytes/50))
end
merge!(metrics,Dict("status"=>"PASS","initial_sha256"=>bytes2hex(sha256(read(initial))),"code_sha256"=>formation_hashes(),"target_limit"=>400.))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true)

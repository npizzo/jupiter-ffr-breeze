include("model.jl")
using Test
source,baseline,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"]
initial_hash=bytes2hex(sha256(read(source)))
@assert initial_hash=="35cb85c5c7f179ac97dce1658f6779cd88e1d39c40b6193292412d250a7a8c6d"
@assert d["t"]==0 && size(u0)==size(v0)==(512,512)
rel(a,b)=sqrt(sum(abs2,a.-b)/max(sum(abs2,b),eps()))
function reaction(q,cf)
    b=cf.b*cf.f .*q.^2;c=-cf.c .*q.^3
    b.-=mean(b);c.-=mean(c);b,c
end
function reference_terms!(m,o,cf)
    u,v=state(m);q=vorticity(o,u,v);qb,qc=reaction(q,cf)
    g=band_curl(o,q,cf.a).+qb.+qc;fu,fv=velocity_from_curl(o,g)
    interior(m.timestepper.Gⁿ.u)[:,:,1].+=fu
    interior(m.timestepper.Gⁿ.v)[:,:,1].+=fv
end
function evolve(u,v;case="biased_plus",steps=30,t=0.,dt=.0025,reference=false,zero_nl=false)
    m,b=setup_polynomial(u,v;case,t,dt,b_override=zero_nl ? 0. : nothing,c_override=zero_nl ? 0. : nothing)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    cb=reference ? (m->reference_terms!(m,b.op,coefficients(case))) : b
    add_callback!(sim,cb,IterationInterval(1);name=:terms,callsite=TendencyCallsite());run!(sim)
    (;m,b,sim)
end
metrics=Dict{String,Any}()
@testset "Polynomial feedback, existing forcing and controlled sign test" begin
    m,b=setup_polynomial(u0,v0;case="biased_plus");o=b.op
    q=vorticity(o,u0,v0)
    @test Float32.(q)==read!(baseline,Array{Float32}(undef,512,512))
    @test !any(o.mask[1,:])
    # Verify full forcing curl against separately evaluated real-space polynomial
    # and the independent original Fourier band implementation.
    errors=Float64[];m3errors=Float64[]
    for case in CASES
        cf=coefficients(case);u,v=case=="biased_minus" ? reflect_state(u0,v0) : (u0,v0)
        mm,bb=setup_polynomial(u,v;case);qq=vorticity(o,u,v);qb,qc=reaction(qq,cf)
        fill!(interior(mm.timestepper.Gⁿ.u),0.);fill!(interior(mm.timestepper.Gⁿ.v),0.);bb(mm)
        fu=copy(dropdims(Array(interior(mm.timestepper.Gⁿ.u));dims=3))
        fv=copy(dropdims(Array(interior(mm.timestepper.Gⁿ.v));dims=3))
        g=vorticity(o,fu,fv);expected=band_curl(o,qq,cf.a).+qb.+qc
        err=rel(g,expected);push!(errors,err);@test err<1e-10
        @test abs(mean(g))<1e-12*maximum(abs,g)
        @test abs(mean(fu))<1e-13 && abs(mean(fv))<1e-13
        direct=3mean(qq.^2 .*qb);identity=3cf.b*cf.f*(mean(qq.^4)-mean(qq.^2)^2)
        push!(m3errors,abs(direct-identity));@test isapprox(direct,identity;rtol=1e-12,atol=1e-12)
        @test cf.f*direct>=0
        @test isapprox(mean(qq.*qc),-cf.c*(mean(qq.^4)-mean(qq)*mean(qq.^3));rtol=1e-12)
        @test mean(qq.*qc)<0
        # Kinetic-energy work from actual momentum tendency matches curl inversion.
        psi=spectral_curl(o,qq,map(x->x==0 ? 0. : -1/x,o.lambda))
        @test isapprox(mean(u.*fu.+v.*fv),-mean(psi.*expected);rtol=1e-11,atol=1e-13)
    end
    ref=evolve(u0,v0;reference=true);fast=evolve(u0,v0)
    ur,vr=state(ref.m);uf,vf=state(fast.m)
    err=max(rel(uf,ur),rel(vf,vr));@test err<1e-10
    um,vm=reflect_state(u0,v0);mir=evolve(um,vm;case="biased_minus")
    um,vm=state(mir.m);ue,ve=reflect_state(uf,vf)
    mirrorerr=max(rel(um,ue),rel(vm,ve));@test mirrorerr<1e-10
    # Turning both new coefficients off must recover the established implementation.
    off=evolve(u0,v0;zero_nl=true);om,ob=setup_fast(u0,v0;case="forced")
    os=Simulation(om;Δt=.0025,stop_iteration=30,verbose=false)
    add_callback!(os,ob,IterationInterval(1);name=:terms,callsite=TendencyCallsite());run!(os)
    ou,ov=state(om);nu,nv=state(off.m);baselineerr=max(rel(nu,ou),rel(nv,ov));@test baselineerr<1e-10
    half=evolve(u0,v0;steps=13);u,v=state(half.m)
    path=joinpath(out,"restart.jld2");jldsave(path;u,v,t=Float64(half.m.clock.time))
    dd=load(path);res=evolve(dd["u"],dd["v"];steps=17,t=dd["t"])
    us,vs=state(res.m);restarterr=max(rel(us,uf),rel(vs,vf));@test restarterr<1e-10
    fine=evolve(u0,v0;steps=60,dt=.00125);uu,vv=state(fine.m)
    dterr=max(rel(uu,uf),rel(vv,vf));@test dterr<1e-3
    fast.sim.stop_iteration=130;GC.gc();perf=@timed run!(fast.sim)
    @test perf.bytes/100<65536 && perf.time/100<.18
    pilots=Dict{String,Any}()
    for case in CASES
        u,v=case=="biased_minus" ? reflect_state(u0,v0) : (u0,v0)
        mm,bb=setup_polynomial(u,v;case)
        sim=Simulation(mm;Δt=.0025,stop_time=1.,verbose=false)
        saved=Float64[];set_stable_step!(sim,bb)
        add_callback!(sim,bb,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
        add_callback!(sim,s->set_stable_step!(s,bb),IterationInterval(10);name=:adaptive)
        add_callback!(sim,s->push!(saved,Float64(s.model.clock.time)),TimeInterval(.25);name=:saved)
        run!(sim)
        @test abs(mm.clock.time-1.)<1e-10
        @test all(t->any(abs.(saved.-t).<1e-8),[.25,.5,.75,1.])
        u,v=state(mm);qend=vorticity(o,u,v)
        cfl=sim.Δt*(maximum(abs,u)+maximum(abs,v))/(2pi/512)
        @test all(isfinite,qend) && abs(mean(qend))<1e-12 && cfl<.5
        pilots[case]=Dict("end_time"=>mm.clock.time,"end_cfl"=>cfl,"dt"=>sim.Δt)
    end
    testdt=stable_step(interior(m.velocities.u),interior(m.velocities.v),512,512,512/2pi,512/2pi,b.b_f,b.c,.1)
    @test testdt<.1 && testdt*(maximum(abs,u0)+maximum(abs,v0))/(2pi/512)<=.3*(1+1e-12)
    merge!(metrics,Dict("curl_max_relative_error"=>maximum(errors),"quadratic_m3_max_absolute_error"=>maximum(m3errors),
        "relative_velocity_error_30_steps"=>err,"relative_mirror_error"=>mirrorerr,"relative_restart_error"=>restarterr,
        "existing_linear_only_relative_error"=>baselineerr,"half_dt_relative_difference"=>dterr,
        "seconds_per_step"=>perf.time/100,"bytes_per_step"=>perf.bytes/100,"pilots"=>pilots,
        "initial_vorticity_rms"=>sqrt(mean(abs2,q)),"identical_initial_saved_vorticity"=>true))
end
merge!(metrics,Dict("status"=>"PASS","a"=>.0002,"b"=>.002,"c"=>.0002,"dt_ceiling"=>.0025,"viscosity"=>.0001,
    "initial_sha256"=>initial_hash,"code_sha256"=>source_hashes()))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)

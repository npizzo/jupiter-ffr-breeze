include("model.jl")
include("reference.jl")
using Test
source,out=ARGS;ispath(out)&&error("Output exists");mkpath(out)
d=load(source);@assert d["t"]==0 && size(d["u"])==(512,512)
initial_hash=bytes2hex(sha256(read(source)))
@assert initial_hash=="35cb85c5c7f179ac97dce1658f6779cd88e1d39c40b6193292412d250a7a8c6d"
old=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
qold=vorticity(old,d["u"],d["v"]);qh=FFTW.fft(qold)
@assert sum(abs2.(qh)[old.k2.>32^2])/sum(abs2,qh)<1e-20
# Corner coordinates coincide on every second point; no downsampling aliases exist.
q=qold[1:2:end,1:2:end]
o=BandOperator(GRID,GRID;kmin=1.,kmax=KCUT,allow_zonal=true)
u0,v0=velocity_from_curl(o,q)
prepared=joinpath(out,"initial.jld2")
jldsave(prepared;u=u0,v=v0,t=0.,initial_sha256=initial_hash,prepared_initial=true)
rel(a,b)=sqrt(sum(abs2,a.-b)/max(sum(abs2,b),eps()))
function actual_force(u,v,c)
    m,b=setup_stress(u,v;case="biased_plus",override=c)
    fill!(interior(m.timestepper.Gⁿ.u),0.);fill!(interior(m.timestepper.Gⁿ.v),0.);b(m)
    fu=copy(dropdims(Array(interior(m.timestepper.Gⁿ.u));dims=3));fv=copy(dropdims(Array(interior(m.timestepper.Gⁿ.v));dims=3))
    fu,fv
end
function evolve(u,v;case="biased_plus",steps=20,t=0.,dt=DTMAX,reference_force=false,coriolis=nothing)
    m,b=setup_stress(u,v;case,t,dt,coriolis)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    function slow!(m)
        uu,vv=state(m);r=reference(uu,vv,coefficients(case))
        interior(m.timestepper.Gⁿ.u)[:,:,1].+=r.fu;interior(m.timestepper.Gⁿ.v)[:,:,1].+=r.fv
    end
    add_callback!(sim,reference_force ? slow! : b,IterationInterval(1);name=:stress,callsite=TendencyCallsite());run!(sim)
    (;m,b,sim)
end
metrics=Dict{String,Any}()
@testset "Conservative projected stress model" begin
    @test rel(vorticity(o,u0,v0),q)<1e-12
    @test isapprox(mean(abs2,q),mean(abs2,qold);rtol=1e-12)
    @test abs(mean(q))<1e-13
    errors=Float64[];workerrors=Float64[]
    for case in CASES
        c=coefficients(case);u,v=case=="biased_minus" ? reflect_state(u0,v0) : (u0,v0)
        fu,fv=actual_force(u,v,c);r=reference(u,v,c)
        err=max(rel(fu,r.fu),rel(fv,r.fv));push!(errors,err);@test err<1e-10
        w=mean(u.*fu.+v.*fv);push!(workerrors,abs(w-r.identity));@test isapprox(w,r.identity;rtol=1e-10,atol=1e-12)
        @test abs(mean(fu))+abs(mean(fv))<1e-12
        @test abs(mean(vorticity(o,fu,fv)))<1e-12
        spec=abs2.(FFTW.fft(fu)).+abs2.(FFTW.fft(fv))
        @test sum(spec[o.k2.>KCUT^2])<1e-22+1e-20sum(spec)
        @test reference(u,v,c;part=:cubic).work<=1e-12
    end
    c=(nu0=-.0002,nu1=0.,nu2=0.,f=100.)
    fu,fv=actual_force(u0,v0,c)
    @test rel(fu,spectral_curl(o,u0,-c.nu0.*o.lambda.*(o.k2.<=KCUT^2)))<1e-10
    @test rel(fv,spectral_curl(o,v0,-c.nu0.*o.lambda.*(o.k2.<=KCUT^2)))<1e-10
    # Exact discrete shear reduction: C = Dyy[nu_T(Z) Z], including both projections.
    h=2pi/GRID;u=[.1sin(3*(j-.5)*h) for i=1:GRID,j=1:GRID];v=zeros(GRID,GRID)
    c=coefficients("biased_plus");fu,fv=actual_force(u,v,c);z=vorticity(o,u,v)
    A=(c.nu0.+c.nu1.*z./c.f.+c.nu2.*(z./c.f).^2).*z
    expected=(circshift(A,(0,-1)).-2A.+circshift(A,(0,1)))./h^2
    @test rel(vorticity(o,fu,fv),expected)<1e-9
    # Independent field with modes beyond the projection, to test the adjoint work identity.
    z=[sin(3i*h+5j*h)+.3cos(37i*h-17j*h) for i=1:GRID,j=1:GRID]
    u,v=velocity_from_curl(o,z);fu,fv=actual_force(u,v,c);r=reference(u,v,c)
    @test isapprox(mean(u.*fu.+v.*fv),r.identity;rtol=1e-10,atol=1e-12)
    fast=evolve(u0,v0);slow=evolve(u0,v0;reference_force=true)
    uf,vf=state(fast.m);us,vs=state(slow.m)
    trajectory=max(rel(uf,us),rel(vf,vs));@test trajectory<1e-10
    um,vm=reflect_state(u0,v0);mir=evolve(um,vm;case="biased_minus")
    um,vm=state(mir.m);ue,ve=reflect_state(uf,vf)
    mirror=max(rel(um,ue),rel(vm,ve));@test mirror<1e-10
    # Constant-f Coriolis has no physical effect on a periodic horizontal 2D flow.
    nof=evolve(u0,v0;coriolis=0.);un,vn=state(nof.m)
    f_error=max(rel(un,uf),rel(vn,vf));@test f_error<1e-10
    first=evolve(u0,v0;steps=7);u,v=state(first.m)
    p=joinpath(out,"restart.jld2");jldsave(p;u,v,t=Float64(first.m.clock.time))
    dd=load(p);second=evolve(dd["u"],dd["v"];steps=13,t=dd["t"]);ur,vr=state(second.m)
    restart=max(rel(ur,uf),rel(vr,vf));@test restart<1e-10
    half=evolve(u0,v0;steps=40,dt=DTMAX/2);uh,vh=state(half.m)
    dt_error=max(rel(uh,uf),rel(vh,vf));@test dt_error<1e-4
    fast.sim.stop_iteration=120;GC.gc();perf=@timed run!(fast.sim)
    @test perf.bytes/100<65536 && perf.time/100<.035
    pilots=Dict{String,Any}()
    for case in CASES
        u,v=case=="biased_minus" ? reflect_state(u0,v0) : (u0,v0)
        m,b=setup_stress(u,v;case);sim=Simulation(m;Δt=DTMAX,stop_time=1.,verbose=false)
        add_callback!(sim,b,IterationInterval(1);name=:stress,callsite=TendencyCallsite())
        add_callback!(sim,s->stable!(s,b),IterationInterval(10);name=:adaptive)
        seen=Float64[];add_callback!(sim,s->push!(seen,m.clock.time),TimeInterval(.25);name=:times)
        run!(sim);u,v=state(m)
        @test all(isfinite,u) && all(isfinite,v) && abs(m.clock.time-1.)<1e-10
        @test all(t->any(abs.(seen.-t).<1e-8),[.25,.5,.75,1.])
        ds=diagnostics(u,v,coefficients(case));@test ds["filtered_rossby_max"]<.5
        pilots[case]=Dict("t"=>m.clock.time,"dt"=>sim.Δt,"max_filtered_Ro"=>ds["filtered_rossby_max"])
    end
    merge!(metrics,Dict("max_force_relative_error"=>maximum(errors),"max_work_absolute_error"=>maximum(workerrors),
      "trajectory_relative_error"=>trajectory,"reflected_f_relative_error"=>mirror,"constant_f_velocity_error"=>f_error,
      "restart_relative_error"=>restart,"half_dt_relative_error"=>dt_error,"seconds_per_step"=>perf.time/100,
      "bytes_per_step"=>perf.bytes/100,"initial_zeta_rms"=>sqrt(mean(abs2,q)),"pilots"=>pilots))
end
merge!(metrics,Dict("status"=>"PASS","initial_sha256"=>initial_hash,"prepared_initial_sha256"=>bytes2hex(sha256(read(prepared))),
 "code_sha256"=>source_hashes(),"compute_location"=>"Unity","grid"=>[GRID,GRID]))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)

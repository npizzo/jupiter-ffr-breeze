include("model.jl")
using Test
source,baseline,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"]
initial_hash=bytes2hex(sha256(read(source)))
@assert initial_hash=="35cb85c5c7f179ac97dce1658f6779cd88e1d39c40b6193292412d250a7a8c6d"
@assert d["t"]==0 && size(u0)==size(v0)==(512,512)
rel(a,b)=sqrt(sum(abs2,a.-b)/max(sum(abs2,b),eps()))
lap(q,dx,dy)=(circshift(q,(1,0)).+circshift(q,(-1,0)).-2q)./dx^2+
            (circshift(q,(0,1)).+circshift(q,(0,-1)).-2q)./dy^2
# Independent, allocating momentum stencil, used only for short validation.
function reference_terms!(m;f=1.)
    u,v=state(m);o=BandOperator(size(u,1),size(u,2);kmin=12.,kmax=16.)
    q=vorticity(o,u,v);qa=sign(f).*min.(sign(f).*q,0.);dx=o.lx/o.nx;dy=o.ly/o.ny
    interior(m.timestepper.Gⁿ.u)[:,:,1] .-= 0.0004 .*(circshift(qa,(0,-1)).-qa)./dy
    interior(m.timestepper.Gⁿ.v)[:,:,1] .+= 0.0004 .*(circshift(qa,(-1,0)).-qa)./dx
end
function evolve(u,v;steps=30,t=0.,f=1.,reference=false,dt=.0025)
    m,b=setup_anticyclone(u,v;t,f,dt)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    cb=reference ? (m->reference_terms!(m;f)) : b
    add_callback!(sim,cb,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
    run!(sim);(;m,b,sim)
end
metrics=Dict{String,Any}()
@testset "Selective diffusion: conservation, dissipation, reflection, restart, performance" begin
    m,b=setup_anticyclone(u0,v0);o=b.op;dx=o.lx/o.nx;dy=o.ly/o.ny
    q=vorticity(o,u0,v0);qa=min.(q,0.)
    @test Float32.(q)==read!(baseline,Array{Float32}(undef,512,512))
    fill!(interior(m.timestepper.Gⁿ.u),0.);fill!(interior(m.timestepper.Gⁿ.v),0.);b(m)
    fu=copy(dropdims(Array(interior(m.timestepper.Gⁿ.u));dims=3))
    fv=copy(dropdims(Array(interior(m.timestepper.Gⁿ.v));dims=3))
    D=vorticity(o,fu,fv);expected=0.0004 .*lap(qa,dx,dy)
    curlerr=rel(D,expected);@test curlerr<1e-12
    @test abs(mean(D))<1e-12*maximum(abs,D)
    @test abs(mean(fu))<1e-14 && abs(mean(fv))<1e-14
    work=mean(u0.*fu.+v0.*fv);exactwork=-.0004*mean(abs2,qa)
    @test work<0 && isapprox(work,exactwork;rtol=1e-12)
    zwork=mean(q.*D)
    fluxwork=-.0004*(mean((circshift(q,(-1,0)).-q).*(circshift(qa,(-1,0)).-qa))/dx^2+
                     mean((circshift(q,(0,-1)).-q).*(circshift(qa,(0,-1)).-qa))/dy^2)
    @test zwork<0 && isapprox(zwork,fluxwork;rtol=1e-12)
    # The selector is zero on cyclonic values, unchanged on anticyclonic ones.
    probe=ComplexF64.([-2. 0.;1. 3.]);select_anticyclone!(probe,1.,1.)
    @test real.(probe)==[-2. 0.;0. 0.]
    probe=ComplexF64.([-2. 0.;1. 3.]);select_anticyclone!(probe,-1.,1.)
    @test real.(probe)==[0. 0.;1. 3.]
    ref=evolve(u0,v0;reference=true);fast=evolve(u0,v0)
    ur,vr=state(ref.m);uf,vf=state(fast.m)
    err=max(rel(uf,ur),rel(vf,vr));@test err<1e-10
    um,vm=reflect_state(u0,v0);mir=evolve(um,vm;f=-1.)
    um,vm=state(mir.m);ue,ve=reflect_state(uf,vf)
    mirrorerr=max(rel(um,ue),rel(vm,ve));@test mirrorerr<1e-10
    half=evolve(u0,v0;steps=13);u,v=state(half.m)
    path=joinpath(out,"restart.jld2");jldsave(path;u,v,t=Float64(half.m.clock.time))
    dd=load(path);res=evolve(dd["u"],dd["v"];steps=17,t=dd["t"])
    us,vs=state(res.m);restarterr=max(rel(us,uf),rel(vs,vf));@test restarterr<1e-10
    fine=evolve(u0,v0;steps=60,dt=.00125);uu,vv=state(fine.m)
    dterr=max(rel(uu,uf),rel(vv,vf));@test dterr<1e-3
    fast.sim.stop_iteration=130;GC.gc();perf=@timed run!(fast.sim)
    @test perf.bytes/100<65536 && perf.time/100<.15
    fast.sim.stop_iteration=800;run!(fast.sim)
    u,v=state(fast.m);qend=vorticity(o,u,v)
    cfl=.0025*(maximum(abs,u)+maximum(abs,v))/dx
    @test isfinite(cfl) && cfl<.5 && abs(fast.m.clock.time-2)<1e-10
    @test .5*(mean(abs2,u)+mean(abs2,v))<.5*(mean(abs2,u0)+mean(abs2,v0))
    @test abs(mean(qend))<1e-12
    merge!(metrics,Dict("curl_relative_error"=>curlerr,"relative_velocity_error_30_steps"=>err,
        "relative_mirror_error"=>mirrorerr,"relative_restart_error"=>restarterr,"half_dt_relative_difference"=>dterr,
        "initial_selective_energy_rate"=>work,"initial_selective_enstrophy_rate"=>zwork,
        "mean_vorticity_tendency"=>mean(D),"seconds_per_step"=>perf.time/100,
        "bytes_per_step"=>perf.bytes/100,"pilot_end_time"=>fast.m.clock.time,"pilot_end_cfl"=>cfl,
        "identical_initial_saved_vorticity"=>true))
end
merge!(metrics,Dict("status"=>"PASS","selective_viscosity"=>.0004,"dt"=>.0025,"viscosity"=>.0001,
    "initial_sha256"=>initial_hash,"code_sha256"=>source_hashes()))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)

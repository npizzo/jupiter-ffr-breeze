include("../viscous2d_long512_v1/model.jl")
include("model.jl")
using Test
source,baseline,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);u0,v0=d["u"],d["v"]
initial_hash=bytes2hex(sha256(read(source)))
@assert initial_hash=="35cb85c5c7f179ac97dce1658f6779cd88e1d39c40b6193292412d250a7a8c6d"
@assert d["t"]==0 && size(u0)==size(v0)==(512,512)
function evolve(u,v;reference=false,steps=30,t=0.)
    if reference
        o=BandOperator(512,512;kmin=24.,kmax=32.,allow_zonal=false)
        m=build_model(o;f=1.);seed!(m,u,v);m.clock.time=t
        b=LinearTerms(o;a=.000175,nu=.0001,dt=.00125)
    else
        m,b=setup_halfscale(u,v;t)
    end
    sim=Simulation(m;Δt=.00125,stop_iteration=steps,verbose=false)
    add_callback!(sim,b,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
    run!(sim);(;m,b,sim)
end
rel(a,b)=sqrt(sum(abs2,a.-b)/sum(abs2,b))
metrics=Dict{String,Any}()
@testset "Matched IC, half wavelength, net growth and efficient startup" begin
    o=BandOperator(512,512;kmin=24.,kmax=32.,allow_zonal=false)
    saved=read!(baseline,Array{Float32}(undef,512,512))
    @test Float32.(vorticity(o,u0,v0))==saved
    # Exact continuum match; compare native discrete rates on paired Fourier modes.
    @test (.000175-.0001)*4 ≈ (.0004-.0001)
    reference=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
    rate_errors=Float64[];ref_rates=Float64[];new_rates=Float64[]
    for j=1:512,i=1:512
        reference.mask[i,j] || continue
        kx=i-1<=256 ? i-1 : i-1-512;ky=j-1<=256 ? j-1 : j-1-512
        ii=mod(2kx,512)+1;jj=mod(2ky,512)+1
        @test o.mask[ii,jj]
        gr=.0004*reference.k2[i,j]-.0001*reference.lambda[i,j]
        gn=.000175*o.k2[ii,jj]-.0001*o.lambda[ii,jj]
        push!(rate_errors,abs(gn-gr)/gr);push!(ref_rates,gr);push!(new_rates,gn)
    end
    @test maximum(rate_errors)<.02
    @test !any(o.mask[1,:])
    metrics["max_relative_discrete_net_growth_difference"]=maximum(rate_errors)
    metrics["reference_discrete_growth_range"]=[minimum(ref_rates),maximum(ref_rates)]
    metrics["halfscale_discrete_growth_range"]=[minimum(new_rates),maximum(new_rates)]
    # Excite the NEW band explicitly: the unchanged initial field occupies 12–16.
    # Check curl of the optimized forcing against the independent spectral expression.
    q=[sin(26*2pi*(i-1)/512+6*2pi*(j-1)/512) for i=1:512,j=1:512]
    ut,vt=velocity_from_curl(o,q);mt,bt=setup_halfscale(ut,vt)
    fill!(interior(mt.timestepper.Gⁿ.u),0.);fill!(interior(mt.timestepper.Gⁿ.v),0.);bt(mt)
    gu=copy(dropdims(Array(interior(mt.timestepper.Gⁿ.u));dims=3));gv=copy(dropdims(Array(interior(mt.timestepper.Gⁿ.v));dims=3))
    expected=band_curl(o,vorticity(o,ut,vt),.000175)
    @test rel(vorticity(o,gu,gv),expected)<1e-10
    ref=evolve(u0,v0;reference=true);fast=evolve(u0,v0)
    ur,vr=state(ref.m);uf,vf=state(fast.m)
    err=max(rel(uf,ur),rel(vf,vr));@test err<1e-10
    half=evolve(u0,v0;steps=13);u,v=state(half.m)
    p=joinpath(out,"restart.jld2");jldsave(p;u,v,t=Float64(half.m.clock.time))
    dd=load(p);res=evolve(dd["u"],dd["v"];steps=17,t=dd["t"])
    us,vs=state(res.m);restart_err=max(rel(us,uf),rel(vs,vf));@test restart_err<1e-10
    fast.sim.stop_iteration=130;GC.gc();perf=@timed run!(fast.sim)
    @test perf.bytes/100<65536
    @test perf.time/100<.15
    fast.sim.stop_iteration=1600;run!(fast.sim)
    u,v=state(fast.m);cfl=.00125*(maximum(abs,u)+maximum(abs,v))/(2pi/512)
    @test isfinite(cfl) && cfl<.5
    @test abs(fast.m.clock.time-2)<1e-10
    merge!(metrics,Dict("relative_velocity_error_30_steps"=>err,"relative_restart_error"=>restart_err,
        "seconds_per_step"=>perf.time/100,"bytes_per_step"=>perf.bytes/100,"pilot_end_time"=>fast.m.clock.time,
        "pilot_end_cfl"=>cfl,"identical_initial_saved_vorticity"=>true,"new_forcing_band_verified"=>true))
end
merge!(metrics,Dict("status"=>"PASS","forcing_coefficient"=>.000175,"dt"=>.00125,"viscosity"=>.0001,
    "initial_sha256"=>initial_hash,"code_sha256"=>source_hashes(),"node"=>get(ENV,"SLURMD_NODENAME","")))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)

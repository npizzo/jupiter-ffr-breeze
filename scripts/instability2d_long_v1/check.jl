include("audit.jl")
include("spectra.jl")
using JLD2, SHA
using Test
out=ARGS[1];ispath(out)&&error("Output exists");mkpath(out)
function evolve(o,u,v;a=0.009375,nu=0.001,alpha=1.,f=1.,dt=0.01,n=20,advection=true,remove_mean_drag=true)
    m=build_model(o;f,advection);seed!(m,u,v)
    b=Audit(o;a,nu,alpha,dt,remove_mean_drag)
    sim=Simulation(m;Δt=dt,stop_iteration=n,verbose=false)
    attach_audit!(sim,b);run!(sim);finish!(b,m)
    state(m),b
end
metrics=Dict{String,Any}()
@testset "Light experiment validation" begin
    o=BandOperator(64,64;kmin=6.,kmax=8.,allow_zonal=false)
    u,v=random_seed(o,1701)
    (up,vp),b=evolve(o,u,v)
    @test b.finished==19
    @test b.max_increment_error<1e-9
    @test b.max_moment_error<1e-9
    @test b.max_drag_work<=0
    @test b.max_viscous_work<=0
    @test b.min_injection_work>=0
    metrics["native_audit"]=report(b)
    (um,vm),_=evolve(o,u,v;f=-1.)
    metrics["f_flip_error"]=max(norm(up-um),norm(vp-vm))/max(norm(up),norm(vp))
    @test metrics["f_flip_error"]<1e-10
    ur,vr=reflect_state(u,v);(ua,va),_=evolve(o,ur,vr;f=-1.)
    ue,ve=reflect_state(up,vp)
    metrics["reflection_error"]=max(norm(ua-ue),norm(va-ve))/max(norm(up),norm(vp))
    @test metrics["reflection_error"]<1e-10
    # Independent analytic cubic-drag decay, uniform velocity, f=0.
    (ud,vd),_=evolve(o,fill(0.2,64,64),zeros(64,64);a=0.,nu=0.,alpha=1.,f=0.,remove_mean_drag=false)
    expected=0.2/sqrt(1+2*0.2^2*0.2)
    metrics["drag_decay_error"]=abs(mean(ud)-expected)
    @test metrics["drag_decay_error"]<1e-9
    # Independent analytic viscosity-only modal decay.
    q=[cos(6*2pi*i/64)*sin(2pi*j/64) for i=0:63,j=0:63]
    uq,vq=velocity_from_curl(o,q)
    (un,vn),_=evolve(o,uq,vq;a=0.,nu=0.001,alpha=0.,f=0.,advection=false)
    expectedq=q.*exp(-0.001*o.lambda[7,2]*0.2)
    metrics["viscous_decay_error"]=norm(vorticity(o,un,vn)-expectedq)/norm(expectedq)
    @test metrics["viscous_decay_error"]<1e-9
end
@testset "Long-pair energy and spectral validation" begin
    source=ARGS[2];d=load(source)
    o=BandOperator(128,128;kmin=6.,kmax=8.,allow_zonal=false)
    u,v=d["u"],d["v"]
    @test size(u)==(128,128)
    @test abs(d["t"]*0.6-120)<1e-8
    for amp in (0.009375,0.)
        (ue,ve),b=evolve(o,u,v;a=amp,n=100)
        @test b.max_energy_error<1e-9
        @test maximum(abs,b.energy)>0
        amp==0 && @test b.energy[2]==0
        u0,v0=velocity_from_curl(o,b.stages[1].q)
        s=spectral_checks(o,b.stages[1].q,b.stages[1].parts,u0,v0)
        @test s.eerr<1e-9
        @test s.terr<1e-9
        @test all(abs.(s.s.T[s.s.k.<6,2]).<1e-12)
        en0=mean(u.^2 .+v.^2)/2;en1=mean(ue.^2 .+ve.^2)/2
        @test abs(en1-en0-sum(b.energy))/en0<1e-9
        metrics[amp==0 ? "unforced_check" : "forced_check"]=Dict("energy_closure"=>b.max_energy_error,"parseval_error"=>s.eerr,"work_error"=>s.terr,"net_advection_ratio"=>s.net)
    end
    # Independent manufactured conservative transfer from k=7 to k=2.
    q=[cos(2*2pi*i/128)+cos(7*2pi*i/128) for i=0:127,j=0:127]
    p=[cos(2*2pi*i/128)-o.lambda[8,1]/o.lambda[3,1]*cos(7*2pi*i/128) for i=0:127,j=0:127]
    z=zeros(128,128);s=spectral_budget(o,q,(p,z,z,z))
    expected=-1/(2o.lambda[3,1])
    @test abs(s.flux[3]-expected)<1e-12
    @test abs(s.flux[8])<1e-12
    metrics["manufactured_inverse_flux"]=s.flux[3]
    metrics["source_checkpoint_sha256"]=bytes2hex(sha256(read(source)))
    metrics["initial_energy"]=mean(u.^2 .+v.^2)/2
end
metrics["status"]="PASS";metrics["compute_location"]="Unity"
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
println("LIGHT_VALIDATION_PASS");TOML.print(stdout,metrics;sorted=true)

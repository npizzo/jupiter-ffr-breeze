include("core.jl")
using .Gamma2026, Test, TOML, JLD2, Statistics, GeophysicalFlows, CUDA
const QG=Gamma2026.QG
mode=ARGS[1];cfg=TOML.parsefile(ARGS[2]);out=ARGS[3];mkpath(out)
small=merge(cfg,Dict("n"=>128,"initial_inverse_wavenumber_m"=>2.0e6,
    "max_dt_seconds"=>30.0,"kind"=>"none","stepper"=>"RK4","annulus_fraction"=>0.1))
function checks(device)
    @testset "analytic wave and invariants" begin
        p=QG.Problem(device;nx=64,Lx=2pi,dt=0.001,stepper="RK4",T=Float64)
        x=Array(p.grid.x);y=Array(p.grid.y)
        q=[-5sin(xi+2yi) for xi in x, yi in y]
        QG.set_q!(p,device_array(device)(q))
        @test QG.energy(p)≈1.25 rtol=1e-12
        e0=QG.energy(p);z0=QG.enstrophy(p);q0=copy(p.sol)
        stepforward!(p,50);QG.updatevars!(p)
        @test maximum(abs,Array(p.sol-q0))<1e-8
        @test QG.energy(p)≈e0 rtol=1e-11
        @test QG.enstrophy(p)≈z0 rtol=1e-11
    end
    @testset "initialization and restart" begin
        p=Gamma2026.problem(small,device);a=Gamma2026.problem(small,device)
        @test QG.energy(p)≈0.5 rtol=1e-12
        @test maximum(abs,Array(p.sol-a.sol))==0
        QG.updatevars!(p);@test abs(mean(Array(p.vars.q)))<1e-11
        for i=1:20;Gamma2026.safe_step!(p,small);end
        ck=joinpath(out,"restart-test.jld2");Gamma2026.checkpoint!(p,small,ck,0.5)
        b=Gamma2026.problem(small,device);Gamma2026.restore!(b,small,ck)
        for i=1:20
            Gamma2026.safe_step!(p,small);Gamma2026.safe_step!(b,small)
        end
        @test maximum(abs,Array(p.sol-b.sol))<1e-10
        @test p.clock.t==b.clock.t
        @test abs(QG.energy(p)/0.5-1)<1e-6
    end
    @testset "polar and flat geometry" begin
        p=Gamma2026.problem(merge(small,Dict("kind"=>"polar")),device)
        f=Gamma2026.problem(merge(small,Dict("kind"=>"flat")),device)
        @test maximum(abs,Array(p.sol-f.sol))==0
        eta=Array(p.params.eta);flat=Array(f.params.eta)
        @test maximum(eta)==0
        @test minimum(flat)≈-cfg["trap_Lgamma"]^2/2
        @test eta != flat
        E0=QG.energy(p)
        for i=1:50;Gamma2026.safe_step!(p,small);end
        @test isfinite(QG.energy(p))
        @test abs(QG.energy(p)/E0-1)<0.02
    end
end
if mode=="unit"
    checks(CPU())
    open(joinpath(out,"unit.toml"),"w") do io
        TOML.print(io,Dict("status"=>"PASSED","source_hash"=>Gamma2026.source_hash()))
    end
elseif mode=="gpu"
    CUDA.functional() || error("GPU not available")
    unit=TOML.parsefile(joinpath(out,"unit.toml"))
    unit["status"]=="PASSED" && unit["source_hash"]==Gamma2026.source_hash() || error("Unit gate missing/stale")
    checks(GPU())
    cpu=Gamma2026.problem(small,CPU());gpu=Gamma2026.problem(small,GPU())
    for i=1:10
        Gamma2026.safe_step!(cpu,small);Gamma2026.safe_step!(gpu,small)
    end
    @test maximum(abs,Array(cpu.sol)-Array(gpu.sol))/maximum(abs,Array(cpu.sol))<1e-8
    p=Gamma2026.problem(cfg,GPU())
    for i=1:20;Gamma2026.safe_step!(p,cfg);end
    CUDA.synchronize();t0=p.clock.t;w0=time()
    for i=1:100;Gamma2026.safe_step!(p,cfg);end
    CUDA.synchronize();elapsed=time()-w0
    seconds= (p.clock.t-t0)*Gamma2026.scales(cfg).T
    estimated=cfg["end_days"]*86400/seconds*elapsed/3600
    m=Gamma2026.metrics(p)
    @test 0.7 <= m.E/0.5 <= 1.02
    estimated<cfg["max_wall_hours"]*0.8 || error("Estimated pilot cost exceeds throughput gate: $estimated hours")
    gate=Dict("status"=>"PASSED","source_hash"=>Gamma2026.source_hash(),"n"=>cfg["n"],
       "measured_steps"=>100,"measured_wall_seconds"=>elapsed,"measured_simulated_seconds"=>seconds,
       "estimated_run_hours"=>estimated,"energy_fraction"=>m.E/0.5,
       "scope"=>"Numerical sanity and throughput only; not resolution convergence or exact paper reproduction")
    open(joinpath(out,"gpu.toml"),"w") do io;TOML.print(io,gate);end
    println("GPU GATE PASSED: estimated pilot hours = ",estimated)
else
    error("mode must be unit or gpu")
end

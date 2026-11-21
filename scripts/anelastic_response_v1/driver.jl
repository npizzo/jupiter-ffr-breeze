function main()
    mode,out=ARGS[1:2];spec=TOML.parsefile("configs/anelastic-response-v1.toml")
    c=merge(TOML.parsefile("configs/anelastic-sequence-v3/baseline.toml"),Dict("initial_dt"=>.001,"max_dt"=>.2,"wall_seconds"=>1500.,"diagnostic_interval"=>1.,"movie_interval"=>5.,"checkpoint_interval"=>100.))
    arch=mode=="preflight" ? CPU() : GPU();arch isa GPU && !CUDA.functional() && error("GPU required")
    isdir(out)&&error("Refuse existing output");mkpath(out)
    if mode=="preflight"
        c=merge(c,Dict("nx"=>16,"ny"=>16,"nz"=>32,"stop_time"=>.2,"diagnostic_interval"=>.1,"movie_interval"=>.1,"checkpoint_interval"=>.2,"max_dt"=>.02,"wall_seconds"=>240.))
        s,_=run_case(c,joinpath(out,"initial");arch)
        cp=only(readdir(joinpath(out,"initial/checkpoints");join=true))
        s=nothing;GC.gc()
        for (name,sgn,n,stride) in [("control",0,16,1),("plus",1,16,1),("minus",-1,16,1),("fine",0,32,1),("no_observer",0,16,0)]
            q=merge(c,Dict("stop_time"=>.4,"nx"=>n,"ny"=>n))
            hook=(sim,p)->restore_and_seed!(sim,p,c,cp;sgn,amplitude=.1,out=joinpath(out,name))
            ss,meta=run_case(q,joinpath(out,name);arch,init_hook=hook,budget_stride=stride)
            @assert meta["final_time"]>=.4-1e-10
            if name=="control";global checked_state=Dict(k=>Array(interior(f)) for (k,f) in pairs(Oceananigans.prognostic_fields(ss.model)));end
            if name=="no_observer"
                @assert all(checked_state[k]==Array(interior(f)) for (k,f) in pairs(Oceananigans.prognostic_fields(ss.model)))
            end
            for path in ["diagnostics.csv";stride>0 ? ["transfer.csv"] : String[]]
                rows=readlines(joinpath(out,name,path));@assert length(rows)>1
                @assert all(all(isfinite,parse.(Float64,split(row,','))) for row in rows[2:end])
            end
            f=load(joinpath(out,name,"frame_00000040.jld2"));@assert size(f["vorticity_slices"],3)==3
            ss=nothing;GC.gc()
        end
        cp2=only(readdir(joinpath(out,"control/checkpoints");join=true))
        restart,_=run_case(merge(c,Dict("stop_time"=>.6)),joinpath(out,"resumed");arch,pickup=cp2,budget_stride=1)
        direct,_=run_case(merge(c,Dict("stop_time"=>.6)),joinpath(out,"direct");arch,pickup=cp,budget_stride=1)
        err=maximum(maximum(abs,Array(interior(f)) .- Array(interior(Oceananigans.prognostic_fields(direct.model)[k]))) for (k,f) in pairs(Oceananigans.prognostic_fields(restart.model)))
        @assert err<1e-8
        write_toml(joinpath(out,"DRIVER_PASS.toml"),Dict("status"=>"PASS","restart_error"=>err,"observer_bitwise_unchanged"=>true,"scope"=>"Actual production driver, both seed signs, conservative refinement, native transfer, output readback, serialized restart"))
        return
    end
    @assert bytes2hex(open(sha256,spec["checkpoint"]))==spec["checkpoint_sha256"]
    mode in ["screen","control","plus","minus"]||error("Unknown mode")
    if mode=="screen"
        for n in [128,256]
            q=merge(c,Dict("nx"=>n,"ny"=>n,"stop_time"=>12020.,"wall_seconds"=>900.))
            name="grid$n";hook=(sim,p)->restore_and_seed!(sim,p,c,spec["checkpoint"];sgn=0,amplitude=.1,out=joinpath(out,name),frame=spec["source_frame"])
            _,meta=run_case(q,joinpath(out,name);arch,init_hook=hook,budget_stride=1)
            @assert meta["final_time"]>=12020-1e-8 "Screen wall limit; no response release"
            GC.gc();CUDA.reclaim()
        end
        write_toml(joinpath(out,"SCREEN_COMPLETE.toml"),Dict("status"=>"COMPLETE_REVIEW_REQUIRED"))
    else
        gate=TOML.parsefile(ARGS[3]);@assert gate["status"]=="PASS"
        q=merge(c,Dict("stop_time"=>14000.));sgn=mode=="plus" ? 1 : mode=="minus" ? -1 : 0
        hook=(sim,p)->restore_and_seed!(sim,p,c,spec["checkpoint"];sgn,amplitude=.1,out=joinpath(out,"parent"),frame=spec["source_frame"])
        _,meta=run_case(q,joinpath(out,"parent");arch,init_hook=hook,budget_stride=100)
        @assert meta["final_time"]>=14000-1e-8 "Partial run; review before extension"
        write_toml(joinpath(out,"COMPLETE.toml"),Dict("status"=>"COMPLETE","final_time"=>meta["final_time"]))
    end
end

include("calibration.jl")
include("budget.jl")
include("remap.jl")

function restore_and_seed!(sim,p,source_config,checkpoint;sgn,amplitude,out,frame=nothing)
    m=sim.model;n=p.nx;coarse=source_config["nx"];arch=Oceananigans.Architectures.architecture(m.grid)
    saved=sim.callbacks[:saved_step].func
    if n==coarse
        set!(sim;checkpoint);sim.Δt=min(saved.dt,.2)
    else
        @assert n==2coarse
        old,_=build(source_config;arch)
        restored=Simulation(old;Δt=.001,stop_time=sim.stop_time)
        ss=SavedStep(.001);add_callback!(restored,ss,IterationInterval(1);name=:saved_step)
        set!(restored;checkpoint);Oceananigans.initialize!(restored)
        before=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(Oceananigans.prognostic_fields(old)))
        for (key,field) in pairs(Oceananigans.prognostic_fields(m))
            a=before[string(key)]
            q=key==:ρu ? prolong_center(prolong_faces(a,1;periodic=true),2;periodic=true) : key==:ρv ? prolong_center(prolong_faces(a,2;periodic=true),1;periodic=true) : prolong_center(prolong_center(a,1;periodic=true,limited=key!=:ρw),2;periodic=true,limited=key!=:ρw)
            back=key==:ρu ? average_pairs(q[1:2:end,:,:],2) : key==:ρv ? average_pairs(q[:,1:2:end,:],1) : average_pairs(average_pairs(q,1),2)
            @assert maximum(abs,back .- a)<1e-12max(maximum(abs,a),1.)
            set!(field,q)
        end
        pressure=Array(interior(old.dynamics.pressure_anomaly))
        set!(m.dynamics.pressure_anomaly,prolong_center(prolong_center(pressure,1;periodic=true),2;periodic=true))
        Oceananigans.restore_prognostic_state!(m.clock,Oceananigans.prognostic_state(old.clock));sim.Δt=min(ss.dt,.1);saved.dt=sim.Δt
        old=nothing;restored=nothing;before=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    end
    Oceananigans.initialize!(sim)
    a=arrays(m)
    if frame!==nothing && n==coarse
        f=load(frame);ub=dropdims(sum(a.rho .* a.u;dims=3) ./ sum(a.rho);dims=3);vb=dropdims(sum(a.rho .* a.v;dims=3) ./ sum(a.rho);dims=3)
        @assert Float32.(ub)==f["ubar"] && Float32.(vb)==f["vbar"]
        @assert maximum(abs,vec(mean(a.theta;dims=(1,2))) .- f["theta_profile"])<1e-12
    end
    dx=p.length/n;dy=p.length/p.ny
    # Corner streamfunction gives exactly divergence-free face increments.
    xx=reshape((0:n-1) .* dx,:,1);yy=reshape((0:p.ny-1) .* dy,1,:)
    psi=-exp.(2.0 .* (cos.(2pi .* (xx .- p.length/2) ./ p.length) .+ cos.(2pi .* (yy .- p.length/2) ./ p.length) .- 2))
    du=-(circshift(psi,(0,-1)) .- psi) ./ dy;dv=(circshift(psi,(-1,0)) .- psi) ./ dx
    zeta=(dv .- circshift(dv,(1,0))) ./ dx .- (du .- circshift(du,(0,1))) ./ dy
    factor=amplitude*abs(p.f)/maximum(zeta);psi .*=factor;du .*=factor;dv .*=factor
    @assert abs(mean(zeta))<1e-12maximum(abs,zeta)
    divseed=(circshift(du,(-1,0)) .- du) ./ dx .+ (circshift(dv,(0,-1)) .- dv) ./ dy
    @assert maximum(abs,divseed)<1e-12
    # Equal isolated seed energy; parent/seed interference is recorded, not hidden.
    cross=mean(dropdims(sum(a.rho .* a.uf;dims=3) ./ sum(a.rho);dims=3) .* du .+ dropdims(sum(a.rho .* a.vf;dims=3) ./ sum(a.rho);dims=3) .* dv)
    if sgn!=0
        set!(m.momentum.ρu,a.mx .+ sgn .* a.rho .* reshape(du,n,p.ny,1))
        set!(m.momentum.ρv,a.my .+ sgn .* a.rho .* reshape(dv,n,p.ny,1))
        pc=.25 .* (psi .+ circshift(psi,(-1,0)) .+ circshift(psi,(0,-1)) .+ circshift(psi,(-1,-1)))
        set!(m.dynamics.pressure_anomaly,Array(interior(m.dynamics.pressure_anomaly)) .+ reshape(sgn*p.f .* pc,n,p.ny,1))
    end
    Oceananigans.TimeSteppers.update_state!(m;compute_tendencies=true)
    a2=arrays(m)
    div=(circshift(a2.mx,(-1,0,0)) .- a2.mx) ./ dx .+ (circshift(a2.my,(0,-1,0)) .- a2.my) ./ dy .+ diff(a2.mz;dims=3) ./ p.dz
    @assert maximum(abs,div)<1e-8
    @assert a.theta==a2.theta && a.wf==a2.wf
    jldsave(joinpath(out,"seed_template.jld2");du,dv,psi,isolated_ke=.5mean(du .^ 2  .+  dv .^ 2),parent_cross_term=cross,sign=sgn)
    write_toml(joinpath(out,"PICKUP_PASS.toml"),Dict("status"=>"PASS","time"=>Float64(m.clock.time),"checkpoint"=>checkpoint,"seed_sign"=>sgn,"seed_peak_zeta_over_f"=>amplitude,"seed_isolated_ke"=>.5mean(du .^ 2  .+  dv .^ 2),"parent_cross_term"=>cross,"mass_divergence"=>maximum(abs,div),"source_n"=>coarse,"target_n"=>n,"balance"=>"Barotropic geostrophic pressure increment; nonlinear projection/adjustment retained; exclude first 50 time units from sustained-response attribution"))
end
main()

function run_branch(cfg,statefile,out;arch=GPU(),pickup=false)
    @assert !ispath(out);mkpath(out)
    s=load(statefile,"state");m,p=build_model(cfg;arch);install_state!(m,s)
    err=state_error(s,extract_state(m));@assert err<1e-11
    sim=Simulation(m;Δt=cfg["max_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cfg["max_dt"])
    function diffusion_limit(sim)
        rmin=minimum(interior(density_field(m)))
        sim.Δt=min(sim.Δt,.08/(max(p.kc/(rmin*CV),1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diffusion_limit,IterationInterval(1);name=:diffusion_limit)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    budget=Budget(m,p);add_callback!(sim,budget,IterationInterval(1);name=:budget)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt;Oceananigans.TimeSteppers.update_state!(m)
    end
    meta=merge(cfg,Dict("status"=>"RUNNING","source_state"=>statefile,"initial_prognostic_error"=>err,
        "energy_note"=>"Incremental residual resets at mapped initial energy; every-step trapezoidal heat and independent one-unit sampling; no energy compensation"))
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io
        println(io,"time,iteration,dt,wall_seconds,energy,kinetic,internal,potential,power,cooling,top_flux,heat_step,heat_sparse,residual_step,residual_sparse,bottom_heat,max_mach,mass_error,min_density,min_temperature")
    end
    started=time()
    function output(sim)
        budget(sim);bs=budget.state;t=Float64(m.clock.time);e=energy(m,p);a=e.a
        P,cooling,top=power(budget,m)
        near(x)=isapprox(x,round(x);atol=1e-7)
        if near((t-bs.start)/cfg["sparse_interval"]) || t>=cfg["stop_time"]-1e-9
            bs=merge(bs,(;sparse_t=t,sparse_p=P,sparse_heat=bs.sparse_heat+.5*(t-bs.sparse_t)*(P+bs.sparse_p)))
            budget.state=bs
        end
        mach=maximum(sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T))
        masserr=abs(mean(a.rho)/bs.mass0-1)
        @assert minimum(a.rho)>0 && minimum(a.T)>0 && masserr<1e-7 && mach<2
        residual=e.E-bs.E0-bs.heat;bottom=16p.flux*(t-bs.start)
        vals=(t,m.clock.iteration,sim.Δt,time()-started,e.E,e.ke,e.ie,e.pe,P,cooling,top,bs.heat,bs.sparse_heat,
            residual,e.E-bs.E0-bs.sparse_heat,bottom,mach,masserr,minimum(a.rho),minimum(a.T))
        @assert all(isfinite,vals)
        open(path,"a") do io;println(io,join(vals,','));end
        if near((t-bs.start)/cfg["movie_interval"]) || t>=cfg["stop_time"]-1e-9
            write_frame(m,p,out,t)
            @printf("t=%.4f dt=%.7g wall=%.1f residual=%.6g bottom=%.5g Mach=%.4g\n",t,sim.Δt,time()-started,residual,bottom,mach);flush(stdout)
        end
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    try
        run!(sim;checkpoint_at_end=true)
        meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    catch e
        meta["status"]="FAILED";meta["reason"]=sprint(showerror,e);rethrow()
    finally
        meta["final_time"]=Float64(m.clock.time);meta["wall_seconds_actual"]=time()-started
        write_toml(joinpath(out,"metadata.toml"),meta)
        write_toml(joinpath(out,"budget_end.toml"),Dict(string(k)=>v for (k,v) in pairs(budget.state)))
    end
    sim,budget
end

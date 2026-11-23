function run_case(cfg,out;arch=GPU(),pickup=false)
    @assert !ispath(out);mkpath(out);m,p=build_model(cfg;arch)
    get(cfg,"probe_amplitude",0.)>0 && add_probe!(m,p,cfg)
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    av=Field(Average(m.velocities.u,dims=1))
    function jetmean(model);compute!(av);interior(p.umean).=interior(av);nothing;end
    jetmean(m);add_callback!(sim,jetmean,IterationInterval(1);name=:jet_mean,callsite=Oceananigans.UpdateStateCallsite())
    cap=min(cfg["max_dt"],.1/max(2abs(p.omega),1e-20),.15p.tau)
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=cap)
    function step_limit(sim)
        nu=maximum(interior(m.closure_fields.νₑ));rmin=minimum(interior(m.dynamics.dry_density))
        diffusivity=max(4nu/3,p.entropy_diffusion_factor*2.5/CV*nu/p.pr,p.kc*max(p.conduction_factor,p.top_factor)/(max(rmin,1e-30)*CV))
        sim.Δt=min(sim.Δt,.08/(max(diffusivity,1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,step_limit,IterationInterval(1);name=:step_limit)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt;Oceananigans.TimeSteppers.update_state!(m);jetmean(m);Oceananigans.TimeSteppers.update_state!(m)
        sim.stop_time=Float64(m.clock.time)+get(cfg,"replay_duration",1.)
    end
    a0=arrays(m);theta0=Array(interior(m.formulation.potential_temperature_density))./a0.rho;tp0=vec(mean(a0.T;dims=(1,2)));mass0=mean(a0.rho)
    limits=guard_spec(cfg);z=[zcell(k,p) for k=1:p.nz];planes=[argmin(abs.(z.-q)) for q in (.2,.5,.9)]
    meta=merge(cfg,Dict("status"=>"RUNNING","start_time"=>Float64(m.clock.time),"guard_spec"=>limits,"native_screen"=>spatial_screen(m,p),"production_release"=>false,"time_unit"=>"nondimensional","initialization"=>p.seed_info,"video_heights"=>z[planes]))
    write_toml(joinpath(out,"metadata.toml"),meta)
    io=get(cfg,"record_stages",false) ? stage_recorder(sim,p,out) : nothing
    function guard(model)
        s=health(model,theta0,tp0,mass0);reason=guard_reason(s,limits)
        if reason!="PASS"
            write_toml(joinpath(out,"GUARD_STOP.toml"),Dict("time"=>Float64(model.clock.time),"reason"=>reason,"health"=>s,"limits"=>limits))
            save_bad(model,out,reason);error("GUARD_STOP: "*reason)
        end
    end
    add_callback!(sim,guard,IterationInterval(1);name=:shared_guard,callsite=Oceananigans.UpdateStateCallsite())
    csvpath=joinpath(out,"diagnostics.csv")
    open(csvpath,"w") do f;println(f,"time,energy,kinetic_energy,external_thermal_power,jet_power,cumulative_external_work,energy_defect,energy_defect_over_initial_ke,Tmin,mach,profile_drift,theta_departure,nonzonal_rms");end
    tr=Ref((t=NaN,E0=NaN,K0=NaN,Q=0.,J=0.,work=0.));started=time()
    function output(sim)
        a=arrays(m);t=Float64(m.clock.time);dv=p.dx*p.dy*p.dz
        rx=.5 .*(a.rho.+circshift(a.rho,(1,0,0)));ry=yface_density(a.rho);rz=.5 .*(a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
        ke=.5dv*(sum(rx.*a.uf.^2)+sum(ry.*a.vf.^2)+sum(rz.*a.wf[:,:,2:end-1].^2))
        E=dv*sum(a.rho.*(CV.*a.T.+GRAV.*reshape(z,1,1,:)))+ke
        target=reshape([jet_target((j-.5)*p.dy-p.ly/2,zcell(k,p),p) for j=1:p.ny,k=1:p.nz],1,p.ny,p.nz)
        meanu=mean(a.uf;dims=1);jet=p.jet_speed==0 ? 0. : dv*sum(rx.*a.uf.*((target.-meanu)./p.jet_tau))
        qtop=p.top_factor*mean(2p.kc.*(a.T[:,:,end].-p.tt)./p.dz)
        cool=p.cap_factor*dv*sum(a.rho.*(a.T.-p.tt).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
        thermal=p.lx*p.ly*(p.bottom_factor*p.flux-qtop)-cool
        old=tr[]
        if !isfinite(old.t);old=(t=t,E0=E,K0=ke,Q=thermal,J=jet,work=0.);end
        work=old.work+.5*(t-old.t)*(old.Q+old.J+thermal+jet);defect=E-old.E0-work
        relative=defect/max(old.K0,1e-8*abs(old.E0))
        tr[]=(t=t,E0=old.E0,K0=old.K0,Q=thermal,J=jet,work=work)
        s=health(m,theta0,tp0,mass0);reason=guard_reason(s,limits);@assert reason=="PASS" reason
        nr=sqrt(mean((a.uf.-meanu).^2)+mean(a.v.^2)+mean(a.w.^2))
        open(csvpath,"a") do f;println(f,join((t,E,ke,thermal,jet,work,defect,relative,s["temperature_min"],s["mach_max"],s["mean_temperature_drift"],s["theta_log_departure"],nr),','));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=sim.stop_time
            dye=Array(interior(m.tracers.dye))./a.rho;zeta=channel_vorticity(a.uf,a.vf,p)
            jldsave(joinpath(out,@sprintf("frame_%07d.jld2",round(Int,100t)));t,z,lx=p.lx,ly=p.ly,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),f=p.f)
            jldsave(joinpath(out,@sprintf("yz_%07d.jld2",round(Int,100t)));t,u=a.uf[1,:,:],v=a.vf[1,:,:],w=a.wf[1,:,:],rho=a.rho[1,:,:],T=a.T[1,:,:],dy=p.dy,dz=p.dz)
        end
        if get(cfg,"save_3d",false) && (t==meta["start_time"]||t>=sim.stop_time)
            jldsave(joinpath(out,@sprintf("full_%07d.jld2",round(Int,100t)));t,z,u_face=a.uf,v_face=a.vf,w_face=a.wf,density=a.rho,temperature=a.T)
        end
        if get(cfg,"enforce_energy",false) && abs(relative)>get(cfg,"energy_limit",.01)
            write_toml(joinpath(out,"ENERGY_STOP.toml"),Dict("time"=>t,"relative_defect"=>relative,"limit"=>get(cfg,"energy_limit",.01)));error("ENERGY_STOP")
        end
        if isapprox(t,round(t);atol=1e-7)||t>=sim.stop_time
            @printf("t=%.4f wall=%.1fs Mach=%.4f Tmin=%.5f profile=%.4g energy/K0=%.4g\n",t,time()-started,s["mach_max"],s["temperature_min"],s["mean_temperature_drift"],relative);flush(stdout)
        end
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    try
        run!(sim;checkpoint_at_end=true)
        meta["status"]=m.clock.time>=sim.stop_time-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    catch err
        save_bad(m,out,sprint(showerror,err));meta["status"]="STOPPED";meta["reason"]=sprint(showerror,err)
        get(cfg,"allow_guard_stop",false) || rethrow()
    finally
        io===nothing || close(io)
        meta["final_time"]=Float64(m.clock.time);meta["wall_seconds_actual"]=time()-started
        write_toml(joinpath(out,"metadata.toml"),meta)
    end
    return sim,meta
end

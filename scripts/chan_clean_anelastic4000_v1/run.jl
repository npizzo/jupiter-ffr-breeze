function run_case(cfg,out;arch=GPU(),pickup=false)
    @assert !ispath(out);mkpath(out)
    m,p=build_model(cfg;arch)
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=cfg["cfl"],max_Δt=min(cfg["max_dt"],.1/(2abs(p.omega))))
    function diffusion_limit(sim)
        rmin=minimum(interior(density_field(m)))
        sim.Δt=min(sim.Δt,.08/(max(p.kc/(rmin*CV),1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diffusion_limit,IterationInterval(1);name=:diffusion_limit)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    function positivity(sim)
        minimum(interior(density_field(m)))>0 && minimum(interior(m.temperature))>0 && minimum(interior(m.formulation.potential_temperature_density))>0 || error("NONPOSITIVE_STATE")
    end
    add_callback!(sim,positivity,IterationInterval(10);name=:positivity)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt
        @assert abs(m.clock.time-cfg["restart_time"])<1e-9
    end
    z=[zcell(k,p) for k=1:p.nz];planes=[argmin(abs.(z.-zz)) for zz in (.2,.5,.9)]
    cz=findall(z .< DC);cap=findall([.5*(z[k]+z[k+1])>.82 for k=1:p.nz-1])
    reference_rho=Array(interior(m.dynamics.reference_state.density));weight=reference_rho./sum(reference_rho;dims=3)
    initial=arrays(m);mass0=mean(initial.rho);Tprofile0=vec(mean(initial.T;dims=(1,2)))
    meta=merge(cfg,Dict("status"=>"RUNNING","initialization"=>p.seed_info,"start_time"=>Float64(m.clock.time),"f"=>p.f,"fh"=>p.fh,"lx"=>4.,"ly"=>4.,"lz"=>1.,"video_heights"=>z[planes],"topology"=>"Periodic, Periodic, Bounded","momentum_forcing"=>"none","rotation_operator"=>"Original periodic MomentumCoriolis / ConstantCartesianCoriolis","cooling"=>"Local -rho*(T-1)*smoothstep(z)/20; no cv multiplier","thermal_equilibrium"=>false,"time_units"=>"H / top isothermal sound speed; nondimensional","energy_note"=>(p.anelastic ? "Compressible energy proxy is not a conserved anelastic energy invariant" : "Staggered KE + internal + gravitational energy minus integrated prescribed heat; unresolved residual is recorded, never compensated")))
    if pickup!==false
        meta["initialization"]=Dict("description"=>"Native checkpoint continuation; no reseeding", "checkpoint"=>pickup)
        meta["diagnostic_baseline"]="Energy proxy ledger and profile_drift reset at continuation start; compare absolute profiles with original run separately"
    end
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io
        println(io,"time,iteration,dt,wall_seconds,kinetic_energy,energy_proxy,external_power,cumulative_heat,energy_defect,defect_over_bottom_heat,max_mach,mass_rms_mach,volume_mach_gt03,volume_mach_gt1,min_density,min_temperature,mass_error,constraint,profile_drift,cap_mean_N2_min,urms_conv,fluctuation_rms_conv,depthmean_ke,zonal_depthmean_ke,tracer_min,tracer_max,cooling_power,top_flux")
    end
    started=time();tracker=Ref((t=NaN,E0=NaN,heat=0.,power=0.,start=Float64(m.clock.time)))
    function fail(reason,t)
        write_toml(joinpath(out,"GUARD_STOP.toml"),Dict("time"=>t,"reason"=>reason));error(reason)
    end
    function output(sim)
        a=arrays(m);t=Float64(m.clock.time);dv=p.dx*p.dy*p.dz
        rx=.5 .*(a.rho.+circshift(a.rho,(1,0,0)));ry=.5 .*(a.rho.+circshift(a.rho,(0,1,0)));rz=.5 .*(a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
        ke=.5dv*(sum(rx.*a.uf.^2)+sum(ry.*a.vf.^2)+sum(rz.*a.wf[:,:,2:end-1].^2))
        E=dv*sum(a.rho.*(CV.*a.T.+GRAV.*reshape(z,1,1,:)))+ke
        qt=mean(2p.kc.*(a.T[:,:,end].-1)./p.dz)
        # Local cooling, including density-temperature covariance; NOT mean-only cooling.
        cooling=dv*sum(a.rho.*(a.T.-1).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
        power=p.heat*(16*(p.flux-qt)-cooling)
        old=tracker[];!isfinite(old.t) && (old=(t=t,E0=E,heat=0.,power=power,start=t))
        heat=old.heat+.5*(t-old.t)*(old.power+power);defect=E-old.E0-heat
        tracker[]=(t=t,E0=old.E0,heat=heat,power=power,start=old.start)
        jldsave(joinpath(out,"budget_sample.jld2");budget=tracker[])
        M=sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T)
        mach=maximum(M);masserr=abs(mean(a.rho)/mass0-1);cons=p.anelastic ? constraint(m,p) : 0.
        drift=maximum(abs,vec(mean(a.T;dims=(1,2)))./Tprofile0 .-1)
        n2=GRAV.*diff(log.(a.theta);dims=3)./p.dz;mean_n2=vec(mean(n2;dims=(1,2)))
        ub=dropdims(sum(a.u.*weight;dims=3);dims=3);vb=dropdims(sum(a.v.*weight;dims=3);dims=3)
        fluctuation=sqrt(mean((a.u[:,:,cz].-ub).^2)+mean((a.v[:,:,cz].-vb).^2)+mean(a.w[:,:,cz].^2))
        urms=sqrt(mean(a.u[:,:,cz].^2 .+a.v[:,:,cz].^2 .+a.w[:,:,cz].^2))
        dye=Array(interior(m.tracers.dye))./a.rho
        vals=(t,m.clock.iteration,sim.Δt,time()-started,ke,E,power,heat,defect,defect/max(16p.flux*(t-old.start),1e-12),mach,sqrt(sum(a.rho.*M.^2)/sum(a.rho)),mean(M.>.3),mean(M.>1.),minimum(a.rho),minimum(a.T),masserr,cons,drift,minimum(mean_n2[cap]),urms,fluctuation,.5mean(ub.^2 .+vb.^2),.5mean(mean(ub;dims=1).^2),minimum(dye),maximum(dye),cooling,qt)
        all(isfinite,vals)||fail("NONFINITE_DIAGNOSTIC",t)
        minimum(a.rho)>0 && minimum(a.T)>0 || fail("NONPOSITIVE_STATE",t)
        masserr<cfg["mass_tolerance"] || fail("MASS_ERROR",t)
        cons<cfg["constraint_tolerance"] || fail("ANELASTIC_CONSTRAINT",t)
        mach<cfg[p.anelastic ? "anelastic_mach_stop" : "compressible_mach_stop"] || fail("MACH_LIMIT",t)
        open(path,"a") do io;println(io,join(vals,','));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            zf=(a.vf.-circshift(a.vf,(1,0,0)))./p.dx.-(a.uf.-circshift(a.uf,(0,1,0)))./p.dy
            zeta=.25 .*(zf.+circshift(zf,(-1,0,0)).+circshift(zf,(0,-1,0)).+circshift(zf,(-1,-1,0)))
            jldsave(joinpath(out,@sprintf("frame_%09d.jld2",round(Int,1000t)));t,z,lx=4.,ly=4.,f=p.f,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),temperature_slices=Float32.(a.T[:,:,planes]),ubar=Float32.(ub),vbar=Float32.(vb),mean_u=Float32.(dropdims(mean(a.u;dims=1);dims=1)),mean_temperature=vec(mean(a.T;dims=(1,2))),mean_density=vec(mean(a.rho;dims=(1,2))),mean_N2=mean_n2,w_rms_profile=vec(sqrt.(mean(a.w.^2;dims=(1,2)))),enthalpy_flux_profile=vec(mean(2.5 .*a.rho.*a.w.*a.T;dims=(1,2))))
        end
        if isapprox(t/cfg["field_interval"],round(t/cfg["field_interval"]);atol=1e-7)||t>=cfg["stop_time"]
            jldsave(joinpath(out,@sprintf("full_%09d.jld2",round(Int,1000t)));t,z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),density=Float32.(a.rho),temperature=Float32.(a.T),theta=Float32.(a.theta))
        end
        @printf("t=%.3f iter=%d dt=%.5g wall=%.1fs Mach=%.4g drift=%.4g capN2=%.4g defect/bottom=%.4g\n",t,m.clock.iteration,sim.Δt,time()-started,mach,drift,minimum(mean_n2[cap]),vals[10]);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(cfg["diagnostic_interval"]);name=:output)
    try
        run!(sim;checkpoint_at_end=true)
        meta["status"]=m.clock.time>=cfg["stop_time"]-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    catch e
        meta["status"]="FAILED";meta["reason"]=sprint(showerror,e);rethrow()
    finally
        meta["final_time"]=Float64(m.clock.time);meta["wall_seconds_actual"]=time()-started;write_toml(joinpath(out,"metadata.toml"),meta)
    end
    sim,meta
end

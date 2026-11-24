function run_case(cfg,out;arch=GPU(),pickup=false)
    @assert !ispath(out);mkpath(out);m,p=build_model(cfg;arch)
    sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["stop_time"],wall_time_limit=cfg["wall_seconds"])
    av=Field(Average(m.velocities.u,dims=1))
    function means(model)
        compute!(av);interior(p.umean).=interior(av)
        interior(p.tmean).=sum(interior(model.temperature);dims=(1,2))./(p.nx*p.ny)
        nothing
    end
    means(m);add_callback!(sim,means,IterationInterval(1);name=:means,callsite=Oceananigans.UpdateStateCallsite())
    conjure_time_step_wizard!(sim;cfl=.25,max_Δt=min(cfg["max_dt"],.1/(2abs(p.omega))))
    function diff_limit(sim)
        rmin=minimum(interior(density_field(m)))
        sim.Δt=min(sim.Δt,.08/(max(p.kc/(rmin*CV),1e-30)*(p.dx^-2+p.dy^-2+p.dz^-2)))
    end
    add_callback!(sim,diff_limit,IterationInterval(1);name=:diff_limit)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(cfg["checkpoint_interval"]),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt;means(m);Oceananigans.TimeSteppers.update_state!(m)
        @assert abs(m.clock.time-cfg["restart_time"])<1e-8
    end
    a0=arrays(m);theta0=Array(interior(m.formulation.potential_temperature));tp0=vec(mean(a0.T;dims=(1,2)));mass0=mean(a0.rho)
    if haskey(cfg,"health_reference_path")
        rr=load(cfg["health_reference_path"]);theta0=rr["theta0"];tp0=rr["tp0"];mass0=rr["mass0"]
    end
    jldsave(joinpath(out,"health_reference.jld2");theta0,tp0,mass0)
    theta0=Oceananigans.Architectures.on_architecture(arch,theta0);tp0=Oceananigans.Architectures.on_architecture(arch,tp0)
    z=[zcell(k,p) for k=1:p.nz];planes=[argmin(abs.(z.-q)) for q in (.2,.5,.9)]
    outside=Oceananigans.Architectures.on_architecture(arch,[!(zcell(k,p)>.9 && abs((j-.5)*p.dy-p.ly/2)>3.5) for i=1:p.nx,j=1:p.ny,k=1:p.nz])
    limits=guard_spec(cfg)
    function check(_)
        model=m
        h=health(model,theta0,tp0,mass0,outside);reason=guard_reason(h,limits)
        if reason!="PASS"
            write_toml(joinpath(out,"GUARD_STOP.toml"),Dict("time"=>Float64(model.clock.time),"reason"=>reason,"health"=>h));error("GUARD_STOP "*reason)
        end
        nothing
    end
    # Full-field guards every 10 outer steps; diagnostics at every saved time.
    # Native finite-state checks remain active. No silent recovery or parameter changes.
    add_callback!(sim,check,IterationInterval(10);name=:guard)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    meta=merge(cfg,Dict("status"=>"RUNNING","start_time"=>Float64(m.clock.time),"initialization"=>p.seed_info,"video_heights"=>z[planes],"time_units"=>"nondimensional","energy_note"=>(p.anelastic ? "Compressible E proxy is not a closed anelastic energy invariant" : "E minus accumulated prescribed thermal and jet work; no compensating heat")))
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,iteration,dt,wall_seconds,kinetic_energy,energy_proxy,external_power,cumulative_work,energy_defect,defect_over_bottom_heat,mach,outside_mach,Tmin,profile_drift,mass_error,constraint,nonzonal_rms,tracer_min,tracer_max,cap_mean_N2_min,jet_rms");end
    spatialpath=joinpath(out,"mach_spatial.csv")
    open(spatialpath,"w") do io
        println(io,"time,mass_rms_mach,volume_gt03,volume_gt05,volume_gt08,volume_gt1,outside_volume_gt05,below09_max_mach,below09_volume_gt05,max_x,max_y,max_z")
    end
    tr=Ref((t=NaN,E0=NaN,work=0.,power=0.,start=Float64(m.clock.time)))
    if haskey(cfg,"budget_reference_path");tr[]=load(cfg["budget_reference_path"],"budget");@assert abs(tr[].t-m.clock.time)<1e-8;end
    started=time()
    function output(sim)
        a=arrays(m);t=Float64(m.clock.time);dv=p.dx*p.dy*p.dz
        rx=.5 .*(a.rho.+circshift(a.rho,(1,0,0)));ry=yface_density(a.rho);rz=.5 .*(a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
        ke=.5dv*(sum(rx.*a.uf.^2)+sum(ry.*a.vf.^2)+sum(rz.*a.wf[:,:,2:end-1].^2))
        E=dv*sum(a.rho.*(CV.*a.T.+GRAV.*reshape(z,1,1,:)))+ke
        meanu=mean(a.uf;dims=1);target=reshape([jet_target((j-.5)*p.dy-p.ly/2,zcell(k,p),p) for j=1:p.ny,k=1:p.nz],1,p.ny,p.nz)
        jet=p.jet_speed==0 ? 0. : dv*sum(rx.*a.uf.*((target.-meanu)./p.jet_tau))
        qt=p.top_factor*mean(2p.kc.*(a.T[:,:,end].-p.tt)./p.dz)
        cool=p.cap_factor*dv*sum(a.rho.*(mean(a.T;dims=(1,2)).-p.tt).*reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
        power=p.lx*p.ly*(p.bottom_factor*p.flux-qt)-cool+jet
        old=tr[];if !isfinite(old.t);old=(t=t,E0=E,work=0.,power=power,start=t);end
        work=old.work+.5*(t-old.t)*(old.power+power);defect=E-old.E0-work
        tr[]=(t=t,E0=old.E0,work=work,power=power,start=old.start)
        jldsave(joinpath(out,"budget_reference.jld2");budget=tr[])
        h=health(m,theta0,tp0,mass0,outside);@assert guard_reason(h,limits)=="PASS"
        cons=p.anelastic ? constraint(m,p) : 0.;@assert !p.anelastic || cons<1e-8
        theta=Array(interior(m.formulation.potential_temperature));n2=GRAV.*diff(log.(theta);dims=3)./p.dz
        cap=findall([.5*(z[k]+z[k+1])>.82 for k=1:p.nz-1]);capn=minimum(vec(mean(n2;dims=(1,2)))[cap])
        M=sqrt.(a.u.^2 .+a.v.^2 .+a.w.^2)./sqrt.((5/3).*a.T)
        outside_host=[!(zcell(k,p)>.9 && abs((j-.5)*p.dy-p.ly/2)>3.5) for i=1:p.nx,j=1:p.ny,k=1:p.nz]
        below=M[:,:,findall(z .<= .9)];imax=Tuple(argmax(M))
        spatial=(t,sqrt(sum(a.rho.*M.^2)/sum(a.rho)),mean(M .> .3),mean(M .> .5),mean(M .> .8),mean(M .> 1.),sum((M .> .5).*outside_host)/sum(outside_host),maximum(below),mean(below .> .5),(imax[1]-.5)*p.dx,(imax[2]-.5)*p.dy-p.ly/2,z[imax[3]])
        @assert all(isfinite,spatial)
        open(spatialpath,"a") do io;println(io,join(spatial,','));end
        dye=Array(interior(m.tracers.dye))./a.rho
        vals=(t,m.clock.iteration,sim.Δt,time()-started,ke,E,power,work,defect,defect/max(p.lx*p.ly*p.flux*(t-old.start),ke,1e-12),h["mach_max"],h["outside_mach"],h["temperature_min"],h["mean_temperature_drift"],h["mass_error"],cons,sqrt(mean((a.uf.-meanu).^2)+mean(a.v.^2)+mean(a.w.^2)),minimum(dye),maximum(dye),capn,sqrt(mean(meanu.^2)))
        @assert all(isfinite,vals)
        open(path,"a") do io;println(io,join(vals,','));end
        if isapprox(t/cfg["movie_interval"],round(t/cfg["movie_interval"]);atol=1e-7)||t>=sim.stop_time
            zeta=channel_vorticity(a.uf,a.vf,p)
            r0=Array(interior(m.dynamics.reference_state.density));weight=r0./sum(r0;dims=3)
            ubar=dropdims(sum(a.u.*weight;dims=3);dims=3);vbar=dropdims(sum(a.v.*weight;dims=3);dims=3)
            jldsave(joinpath(out,@sprintf("frame_%09d.jld2",round(Int,1000t)));t,z,lx=p.lx,ly=p.ly,f=p.f,video_heights=z[planes],vorticity_slices=Float32.(zeta[:,:,planes]),tracer_slices=Float32.(dye[:,:,planes]),ubar=Float32.(ubar),vbar=Float32.(vbar),mean_u=Float32.(dropdims(meanu;dims=1)),mean_temperature=vec(mean(a.T;dims=(1,2))),mean_density=vec(mean(a.rho;dims=(1,2))),mean_N2=vec(mean(n2;dims=(1,2))))
        end
        if cfg["save_3d"] && (isapprox(t/cfg["field_interval"],round(t/cfg["field_interval"]);atol=1e-7)||t>=sim.stop_time)
            jldsave(joinpath(out,@sprintf("full_%09d.jld2",round(Int,1000t)));t,z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),density=Float32.(a.rho),temperature=Float32.(a.T),theta=Float32.(theta))
        end
        @printf("t=%.4f step=%d dt=%.5g wall=%.1fs Mach=%.3f outside=%.3f Tdrift=%.3g defect/bottom=%.3g\n",t,m.clock.iteration,sim.Δt,time()-started,h["mach_max"],h["outside_mach"],h["mean_temperature_drift"],vals[10]);flush(stdout)
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

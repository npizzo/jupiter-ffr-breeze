function run_case(c,ratio,out;arch=GPU(),pickup=false)
    @assert !ispath(out);mkpath(out)
    m,p,info=build(c,ratio;arch)
    τ=p.seed.turnover;tstop=c["turnovers"]*τ
    maxdt=min(.1/p.f,.15/sqrt(p.N2),.08/(p.seed.ν*(2/p.dx^2+1/p.dz^2)))
    get(c,"test",false) && (maxdt=c["initial_dt"])
    sim=Simulation(m;Δt=c["initial_dt"],stop_time=tstop,wall_time_limit=c["wall_seconds"])
    conjure_time_step_wizard!(sim;cfl=c["cfl"],max_Δt=maxdt)
    saved=SavedStep(sim.Δt);add_callback!(sim,saved,IterationInterval(1);name=:saved_step)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(c["checkpoint_turnovers"]*τ),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    if pickup!==false
        set!(sim;checkpoint=pickup);sim.Δt=saved.dt
        Oceananigans.TimeSteppers.update_state!(m)
    end
    meta=merge(c,Dict("status"=>"RUNNING","initialization"=>info,"start_time"=>Float64(m.clock.time),"stop_time"=>tstop,"equations"=>"Breeze dry anelastic; vertical constant rotation; no beta, no jet, no heating or cooling","boundary_conditions"=>"Periodic x/y; impermeable stress-free top/bottom; no scalar flux","initial_balance"=>"Common first baroclinic velocity shape, geostrophic plus nonlinear pressure correction, exact dry buoyancy inversion; finite discrete balance error reported","scope"=>"Internal-deformation evolution of initially mixed-sign eddies. Not convection-driven formation and not a test of compressible Lamb modes.","length_convention"=>"L = sqrt(horizontal face velocity variance / native vorticity variance); U is horizontal face RMS at max vertical-mode amplitude","energy_note"=>"KE and small-amplitude APE proxy are diagnostics, not an exact anelastic invariant","time_units"=>"Nondimensional; turnover = L/U", "video_heights"=>p.z[[argmin(abs.(p.z.-z)) for z in (.25,.6,.9)]]))
    write_toml(joinpath(out,"metadata.toml"),meta)
    path=joinpath(out,"diagnostics.csv")
    open(path,"w") do io;println(io,"time,turnovers,iteration,dt,wall_seconds,kinetic_energy,ape_proxy,max_mach,min_temperature,constraint,theta_mass_relative_error,tracer_mass_relative_error,min_N2,mean_min_N2,w_rms,urms,upper_skew,upper_anti_enstrophy,wall_w_max");end
    started=time();a0=fields_cpu(m);θmass0=sum(a0.rho.*a0.θ);cmass0=sum(a0.rho.*a0.dye)
    planes=[argmin(abs.(p.z.-z)) for z in (.25,.6,.9)]
    function output(sim)
        a=fields_cpu(m);t=Float64(m.clock.time);s=t/τ;dv=p.dx^2*p.dz
        ζf=(a.vf-circshift(a.vf,(1,0,0)))/p.dx-(a.uf-circshift(a.uf,(0,1,0)))/p.dx
        ζ=0.25 .*(ζf+circshift(ζf,(-1,0,0))+circshift(ζf,(0,-1,0))+circshift(ζf,(-1,-1,0)))
        q=ζ[:,:,planes[end]]/p.f
        b=p.g.*(1 .- reshape(p.θ0,1,1,:)./a.θ)
        ke=.5dv*sum(a.rho.*(a.u.^2+a.v.^2+a.w.^2));ape=.5dv*sum(a.rho.*b.^2)/p.N2
        mach=maximum(sqrt.(a.u.^2+a.v.^2+a.w.^2)./sqrt.((5/3).*a.T))
        cons=constraint(m,p);n2=p.g.*diff(log.(a.θ);dims=3)/p.dz
        te=abs(sum(a.rho.*a.θ)/θmass0-1);ce=abs(sum(a.rho.*a.dye)/cmass0-1)
        skew=mean(q.^3)/max(mean(q.^2)^1.5,1e-30);share=sum(q[q.<0].^2)/max(sum(q.^2),1e-30)
        vals=(t,s,m.clock.iteration,sim.Δt,time()-started,ke,ape,mach,minimum(a.T),cons,te,ce,minimum(n2),minimum(mean(n2;dims=(1,2))),sqrt(mean(a.w.^2)),sqrt(mean(a.u.^2+a.v.^2)),skew,share,maximum(abs,a.wf[:,:,[1,end]]))
        reason= !all(isfinite,vals) ? "NONFINITE" : minimum(a.T)<=0 ? "NONPOSITIVE_T" : cons>c["constraint_tolerance"] ? "CONSTRAINT" : mach>c["mach_stop"] ? "MACH_LIMIT" : max(te,ce)>1e-7 ? "SCALAR_CONSERVATION" : nothing
        if reason!==nothing
            write_toml(joinpath(out,"GUARD_STOP.toml"),Dict("reason"=>reason,"time"=>t));error(reason)
        end
        open(path,"a") do io;println(io,join(vals,','));end
        due(interval)=isapprox(s/interval,round(s/interval);atol=1e-7)||t>=tstop-1e-8
        if due(c["movie_turnovers"])
            jldsave(joinpath(out,@sprintf("frame_%09d.jld2",round(Int,1000t)));t,turnovers=s,z=p.z,lx=c["lx"],ly=c["lx"],f=p.f,video_heights=p.z[planes],vorticity_slices=Float32.(ζ[:,:,planes]),tracer_slices=Float32.(a.dye[:,:,planes]),vertical_velocity_slices=Float32.(a.w[:,:,planes]),buoyancy_slices=Float32.(b[:,:,planes]),mean_N2=vec(mean(n2;dims=(1,2))),mean_temperature=vec(mean(a.T;dims=(1,2))))
        end
        if due(c["full_turnovers"])
            jldsave(joinpath(out,@sprintf("full_%09d.jld2",round(Int,1000t)));t,turnovers=s,z=p.z,u_face=Float32.(a.uf),v_face=Float32.(a.vf),w_face=Float32.(a.wf),theta=Float32.(a.θ),temperature=Float32.(a.T),density=Float32.(a.rho),reference_theta=p.θ0,tracer=Float32.(a.dye))
        end
        @printf("t=%.3f turnovers=%.3f iter=%d dt=%.5g wall=%.1fs Mach=%.4g constraint=%.3g minN2=%.4g anti_upper=%.3f\n",t,s,m.clock.iteration,sim.Δt,time()-started,mach,cons,minimum(n2),share);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(c["diagnostic_turnovers"]*τ);name=:output)
    try
        run!(sim;checkpoint_at_end=true)
        meta["status"]=m.clock.time>=tstop-1e-8 ? "COMPLETE" : "WALL_LIMIT"
    catch e
        meta["status"]="FAILED";meta["reason"]=sprint(showerror,e);rethrow()
    finally
        meta["final_time"]=Float64(m.clock.time);meta["final_turnovers"]=m.clock.time/τ;meta["wall_seconds_actual"]=time()-started;write_toml(joinpath(out,"metadata.toml"),meta)
    end
    sim,meta
end

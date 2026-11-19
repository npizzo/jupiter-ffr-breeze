using Breeze, Oceananigans, CUDA, JLD2, TOML, Statistics, Printf, SHA
using Oceananigans: interior
import Breeze.AtmosphereModels as AM
include("budget.jl")
mutable struct PlumeTimeStep
    value::Float64
end
(s::PlumeTimeStep)(sim)=(s.value=sim.Δt;nothing)
Oceananigans.prognostic_state(s::PlumeTimeStep)=(;value=s.value)
Oceananigans.restore_prognostic_state!(s::PlumeTimeStep,x::NamedTuple)=(s.value=x.value;s)

function make_plume(cfg,arch)
    nx,ny,nz=cfg["nx"],cfg["ny"],cfg["nz"]
    L=cfg["lx"];H=40000.;z0=70000.;zi=10000.;g=25.;cp=13000.;θ0=170.
    grid=RectilinearGrid(arch,Float64;size=(nx,ny,nz),x=(-L/2,L/2),y=(-L/2,L/2),z=(0,H),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    R=8.314462618/.0023
    π70=(2e6/1e5)^(R/cp)-g*z0/(cp*θ0)
    @assert π70>0
    p70=1e5*π70^(cp/R)
    theta_ref(z)=θ0*exp(cfg["cap_n2"]*max(z-zi,0)/g)
    constants=ThermodynamicConstants(Float64;gravitational_acceleration=g,dry_air_molar_mass=.0023,dry_air_heat_capacity=cp)
    ref=ReferenceState(grid,constants;base_pressure=p70,potential_temperature=theta_ref)
    f=2*.00017585*sind(55.)
    # Same diffusivity profile as the large-box case, evaluated at the cropped height.
    function kappa(x,y,z,t)
        s=clamp((z + 70000.0 - 60000.0)/20000.0,0.0,1.0)
        1000.0 + (1.0 - 1000.0)*s*s*(3-2s)
    end
    model=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref),thermodynamic_constants=constants,coriolis=FPlane(f=f),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=(ρθ=WENO(order=5),ρqᵛ=WENO(order=5,bounds=(0,1)),dye=WENO(order=5,bounds=(0,1))),closure=ScalarDiffusivity(ν=1000.,κ=(ρθ=kappa,ρqᵛ=kappa,dye=100.)))
    x=((1:nx).-0.5).*L/nx.-L/2;y=((1:ny).-0.5).*L/ny.-L/2;z=((1:nz).-0.5).*H/nz
    pulse(x,y,z)=hypot(x,y)<40000. && abs(z-5000.)<5000. ? cospi(hypot(x,y)/80000.)^2*cospi((z-5000.)/10000.)^2 : 0.
    θ=[theta_ref(zz)+cfg["pulse_kelvin"]*pulse(xx,yy,zz) for xx in x,yy in y,zz in z]
    set!(model;θ,u=0,v=0,w=0)
    rho=vec(Array(interior(ref.density)))
    C=[pulse(xx,yy,zz) for xx in x,yy in y,zz in z]
    set!(model.tracers.dye,C.*reshape(rho,1,1,:))
    (;model,rho,z=z.+z0,x,y,f,p70,theta_ref=theta_ref.(z))
end

function run_plume(cfg,out;arch=GPU(),stop_override=nothing,pickup=false)
    isdir(out)&&!isempty(readdir(out))&&error("Refusing nonempty output directory")
    mkpath(out);mkpath(joinpath(out,"frames"));mkpath(joinpath(out,"fields"))
    arch isa GPU && !CUDA.functional() && error("CUDA unavailable; no CPU fallback")
    q=make_plume(cfg,arch);m=q.model
    stop=isnothing(stop_override) ? cfg["stop_seconds"] : stop_override
    sim=Simulation(m;Δt=cfg["dt"],stop_time=stop,wall_time_limit=cfg["wall_time_limit"])
    conjure_time_step_wizard!(sim;cfl=.2,max_Δt=cfg["dt"])
    saved_dt=PlumeTimeStep(sim.Δt)
    add_callback!(sim,saved_dt,IterationInterval(1);name=:saved_dt)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    b=PlumeBudget(m,out,q.rho,q.z,cfg["lx"]/cfg["nx"])
    add_callback!(sim,b,IterationInterval(1);name=:layer_vorticity_budget,callsite=TendencyCallsite())
    tf=Field(liquid_ice_potential_temperature(m));dz=40000/cfg["nz"];h=cfg["lx"]/cfg["nx"]
    ids=[argmin(abs.(q.z.-z)) for z in (75000.,80500.,82000.)]
    metadata=merge(cfg,Dict("status"=>"RUNNING","f"=>q.f,"reference_base_pressure"=>q.p70,"video_heights_m"=>q.z[ids],"script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"budget_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"budget.jl")))),"gpu"=>arch isa GPU ? CUDA.name(CUDA.device()) : "CPU","model"=>"dry anelastic; periodic horizontal; free-slip impermeable vertical boundaries; warm initial pulse only; no thermostat/sponge/bottom heat/shear/beta","budget"=>"every SSPRK3 stage; density-weighted layer-mean vertical vorticity; native horizontal and vertical advection, Coriolis, viscosity; pressure curl checked; no horizontal forcing","breeze_pin"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b","neutral_control"=>"same bottom pressure and pulse, rebuilt neutral hydrostatic reference","interpretation"=>"isolated plume production test, not LSV formation or a converged solution"))
    writemeta()=open(io->TOML.print(io,metadata),joinpath(out,"metadata.toml"),"w")
    writemeta();started=time();lastsave=Ref(-Inf)
    open(joinpath(out,"diagnostics.csv"),"w") do io
        println(io,"time,wall_seconds,dt,w_max,u_rms,plume_Ro,w_over_N_cells,tracer_top_fraction,tracer_side_fraction,max_layer_step_residual")
    end
    function output(sim)
        t=Float64(m.clock.time);t<=lastsave[]+1e-9&&return;lastsave[]=t
        save_budget!(b,m);compute!(tf)
        uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
        u=(uf.+circshift(uf,(-1,0,0)))./2;v=(vf.+circshift(vf,(0,-1,0)))./2;w=(wf[:,:,1:end-1].+wf[:,:,2:end])./2
        theta=Array(interior(tf));C=Array(interior(m.tracers.dye))./reshape(q.rho,1,1,:)
        ζ=(circshift(v,(-1,0,0)).-circshift(v,(1,0,0)))./(2h).-(circshift(u,(0,-1,0)).-circshift(u,(0,1,0)))./(2h)
        div=(circshift(uf,(-1,0,0)).-uf)./h.+(circshift(vf,(0,-1,0)).-vf)./h
        nearcap=findall(78000. .<= q.z .<= 82000.);wr=maximum(abs,w[:,:,nearcap])
        ur=sqrt(mean(u.^2 .+ v.^2));Ro=ur/(q.f*40000.)
        massC=C.*reshape(q.rho,1,1,:);total=sum(massC)
        top=sum(massC[:,:,q.z.>105000.])/max(total,eps())
        side=repeat(reshape((abs.(q.x).>80000.),:,1),1,length(q.y)).|repeat(reshape((abs.(q.y).>80000.),1,:),length(q.x),1)
        sidefraction=sum(sum(massC;dims=3)[:,:,1][side])/max(total,eps())
        jldsave(joinpath(out,"frames",@sprintf("frame_%06d.jld2",round(Int,t)));t,z=q.z[ids],f=q.f,x=q.x,y=q.y,vorticity=Float32.(ζ[:,:,ids]),tracer=Float32.(C[:,:,ids]),w=Float32.(w[:,:,ids]),divergence=Float32.(div[:,:,ids]))
        if t==0 || abs(t/600-round(t/600))<1e-7 || t>=stop
            jldsave(joinpath(out,"fields",@sprintf("field_%06d.jld2",round(Int,t)));t,z=q.z,x=q.x,y=q.y,rho=q.rho,u_face=Float32.(uf),v_face=Float32.(vf),w_face=Float32.(wf),theta=Float32.(theta),tracer=Float32.(C))
        end
        # Fast, every-minute vertical profiles retain the interface event.
        jldsave(joinpath(out,"frames",@sprintf("profiles_%06d.jld2",round(Int,t)));t,z=q.z,theta_mean=vec(mean(theta;dims=(1,2))),w_rms=vec(sqrt.(mean(w.^2;dims=(1,2)))),zeta_rms=vec(sqrt.(mean(ζ.^2;dims=(1,2)))),theta_reference=q.theta_ref)
        N=sqrt(cfg["cap_n2"]);resolution=N>0 ? wr/N/dz : -1.
        open(joinpath(out,"diagnostics.csv"),"a") do io;println(io,join((t,time()-started,sim.Δt,maximum(abs,w),ur,Ro,resolution,top,sidefraction,b.max_residual),','));end
        @printf("PLUME_PROGRESS case=%s t=%.1f/%.1f wall=%.1f wmax=%.4g Ro=%.4g penetration_cells=%.2f top_tracer=%.3g budget=%.3g\n",cfg["name"],t,stop,time()-started,maximum(abs,w),Ro,resolution,top,b.max_residual);flush(stdout)
        all(isfinite,theta)&&all(isfinite,w)||error("Nonfinite state")
    end
    add_callback!(sim,output,TimeInterval(cfg["output_seconds"]);name=:publication_and_audit)
    sim.output_writers[:checkpointer]=Checkpointer(m;schedule=TimeInterval(1800.),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=true)
    if pickup!==false;set!(sim;checkpoint=pickup);sim.Δt=saved_dt.value;end
    Oceananigans.TimeSteppers.update_state!(m;compute_tendencies=true)
    metadata["segment_start"]=Float64(m.clock.time)
    output(sim);run!(sim;checkpoint_at_end=true);output(sim)
    metadata["status"]=m.clock.time>=stop-1e-6 ? "COMPLETE" : "PARTIAL_WALL_LIMIT"
    metadata["final_time"]=Float64(m.clock.time);metadata["budget_steps"]=b.checks;metadata["max_step_absolute_residual"]=b.max_residual
    writemeta();return (;sim,b,q)
end

if abspath(PROGRAM_FILE)==@__FILE__
    cfg=TOML.parsefile(ARGS[1]);out=ARGS[2];run_plume(cfg,out;arch=length(ARGS)>=3&&ARGS[3]=="CPU" ? CPU() : GPU())
end

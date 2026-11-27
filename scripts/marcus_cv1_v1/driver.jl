include("model.jl")
function sample(m,c,p,info,out)
    uf=Array(interior(m.velocities.u));vf=Array(interior(m.velocities.v));wf=Array(interior(m.velocities.w))
    u=.5*(uf+circshift(uf,(-1,0,0)));v=.5*(vf+circshift(vf,(0,-1,0)));w=.5*(wf[:,:,1:end-1]+wf[:,:,2:end])
    ρ=reshape(info.rho,1,1,:);θ=Array(interior(m.formulation.potential_temperature));T=Array(interior(m.temperature))
    dye=Array(interior(m.tracers.dye))./ρ
    dx=c["lx"]/c["nx"];dy=c["ly"]/c["ny"];dz=(c["zmax"]-c["zmin"])/c["nz"]
    ζ=(vf-circshift(vf,(1,0,0)))/dx-(uf-circshift(uf,(0,1,0)))/dy
    ζ=.25*(ζ+circshift(ζ,(-1,0,0))+circshift(ζ,(0,-1,0))+circshift(ζ,(-1,-1,0)))
    Mx=Array(interior(m.momentum.ρu));My=Array(interior(m.momentum.ρv));Mz=Array(interior(m.momentum.ρw))
    div=(circshift(Mx,(-1,0,0))-Mx)/dx+(circshift(My,(0,-1,0))-My)/dy+diff(Mz;dims=3)/dz
    den=max(maximum(abs,Mx)/dx,maximum(abs,My)/dy,maximum(abs,Mz)/dz,1e-12)
    cons=maximum(abs,div)/den
    rate=Array(interior(p.rate));Ub=Array(interior(p.ubg));Thb=Array(interior(p.thbg));Pi=(reshape(vec(Array(interior(m.dynamics.reference_state.pressure))),1,1,:)/c["p0"]).^(c["R"]/c["cp"])
    # Cell-centred diagnostic energy; convergence needed before any exact-budget claim.
    E=sum(ρ.*(.5*(u.^2+v.^2+w.^2)+c["cp"].*Pi.*θ))*dx*dy*dz
    sink=sum(-rate.*ρ.*(u.*(u-Ub)+v.^2+w.^2+c["cp"].*Pi.*(θ-Thb)))*dx*dy*dz
    vals=(Float64(m.clock.time),maximum(abs,w),maximum(abs,u),cons,minimum(T),E,sink)
    @assert all(isfinite,vals) && minimum(T)>0 && cons<1e-8
    open(joinpath(out,"diagnostics.csv"),"a") do io;println(io,join(vals,','));end
    k=[max(1,round(Int,c["nz"]*q)) for q in (.25,.5,.75)]
    jldsave(joinpath(out,@sprintf("frame_%09d.jld2",round(Int,1000m.clock.time)));t=Float64(m.clock.time),z=info.z,video_heights=info.z[k],vorticity_slices=ζ[:,:,k],tracer_slices=dye[:,:,k],temperature_slices=T[:,:,k],vertical_velocity_slices=w[:,:,k],theta_slices=θ[:,:,k],u_face=uf,v_face=vf,w_face=wf)
    println("FIXTURE t=",m.clock.time," max_w=",vals[2]," constraint=",cons);flush(stdout)
    vals
end
function run_case(c,out;arch=CPU(),vortex=true,rest=false,pickup=false)
    @assert !ispath(out);mkpath(out);m,p,info=build(c;arch,vortex,rest)
    sim=Simulation(m;Δt=c["dt"],stop_time=c["stop_time"],wall_time_limit=300.)
    sim.output_writers[:checkpoint]=Checkpointer(m;schedule=TimeInterval(1.),dir=joinpath(out,"checkpoints"),cleanup=true,verbose=false)
    if pickup!==false;set!(sim;checkpoint=pickup);end
    open(joinpath(out,"diagnostics.csv"),"w") do io;println(io,"time,max_w,max_u,mass_flux_divergence_relative,min_temperature,diagnostic_energy,sponge_power");end
    add_callback!(sim,s->sample(s.model,c,p,info,out),TimeInterval(.5);name=:output)
    run!(sim;checkpoint_at_end=true)
    write_toml(joinpath(out,"metadata.toml"),Dict("status"=>"COMPLETE","fixture_only"=>true,"is_CV1_reproduction"=>false,"final_time"=>m.clock.time,"grid"=>[c["nx"],c["ny"],c["nz"]],"parameters"=>c))
    sim
end
function readback(out)
    rows=readlines(joinpath(out,"diagnostics.csv"));@assert length(rows)>1
    @assert all(all(isfinite,parse.(Float64,split(s,','))) for s in rows[2:end])
    for f in filter(x->startswith(x,"frame_"),readdir(out;join=false))
        a=load(joinpath(out,f));@assert size(a["vorticity_slices"],3)==3
        @assert maximum(abs,a["w_face"][:,:,[1,end]])==0
        for key in ("vorticity_slices","tracer_slices","temperature_slices","vertical_velocity_slices")
            @assert all(isfinite,a[key])
        end
    end
end
function native_buoyancy_check(c)
    m,p,info=build(c;arch=CPU(),rest=true);errs=Float64[];diffs=Float64[]
    for δ in (-.05,-.01,.01,.05)
        set!(m;θ=reshape(info.θref,1,1,:).*ones(c["nx"],c["ny"],1).*(1+δ),qᵗ=0.)
        Oceananigans.TimeSteppers.update_state!(m)
        F=Oceananigans.fields(m);k=c["nz"]÷2
        bn=Breeze.AtmosphereModels.buoyancy_forceᶜᶜᶜ(1,1,k,m.grid,m.dynamics,m.temperature,zeros(c["nx"],c["ny"],c["nz"]),m.microphysics,m.microphysical_fields,m.thermodynamic_constants)
        bc=buoyancy_correction(1,1,k,m.grid,F,p);target=c["g"]*info.rho[k]*δ
        push!(errs,abs(bn+bc-target)/abs(target));push!(diffs,abs(bn-target)/abs(target))
    end
    @assert maximum(errs)<1e-10
    Dict("corrected_buoyancy_max_relative_error"=>maximum(errs),"native_buoyancy_max_relative_difference"=>maximum(diffs),"anomaly_samples"=>[-.05,-.01,.01,.05])
end
function main()
    mode,out=ARGS;@assert mode in ("preflight","gpu_check");@assert !ispath(out);mkpath(out)
    c=fixture();arch=mode=="preflight" ? CPU() : GPU();result=Dict{String,Any}("fixture_only"=>true,"CV1_production_released"=>false)
    mode=="preflight" && merge!(result,native_buoyancy_check(c))
    for (name,vortex,rest) in (("rest",false,true),("shear",false,false),("vortex",true,false))
        a=run_case(c,joinpath(out,name);arch,vortex,rest);readback(joinpath(out,name))
        if rest
            speed=maximum(maximum(abs,Array(interior(f))) for f in values(a.model.velocities));@assert speed<1e-8;result["rest_max_speed"]=speed
        end
        if vortex
            first=run_case(merge(c,Dict("stop_time"=>1.)),joinpath(out,"first");arch)
            cp=only(readdir(joinpath(out,"first/checkpoints");join=true))
            b=run_case(c,joinpath(out,"restart");arch,pickup=cp)
            fa=Oceananigans.prognostic_fields(a.model);fb=Oceananigans.prognostic_fields(b.model)
            err=maximum(maximum(abs,Array(interior(fa[k]))-Array(interior(fb[k])))/max(maximum(abs,Array(interior(fa[k]))),1.) for k in keys(fa))
            @assert err<1e-10;result["restart_relative_error"]=err;readback(joinpath(out,"first"));readback(joinpath(out,"restart"))
        end
        a=nothing;GC.gc();arch isa GPU && CUDA.reclaim()
    end
    result["status"]="PASS";mode=="gpu_check" && (result["gpu"]=CUDA.name(CUDA.device()))
    write_toml(joinpath(out,"DRIVER_PASS.toml"),result);println("MARCUS_FIXTURE_PASS");flush(stdout)
end
main()

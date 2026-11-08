using Breeze, Oceananigans, CUDA, JLD2, TOML, Statistics, Printf, Test, SHA
using Oceananigans: interior
include("validate_jovian_moist_thermodynamics.jl")
root=length(ARGS)>0 ? ARGS[1] : "data/moist-validation-v1"
isdir(root) && !isempty(readdir(root)) && error("Nonempty output root")
mkpath(root)
validate_jovian_moist_thermodynamics(joinpath(root,"thermodynamics"))
CUDA.functional() || error("CUDA required")

function plume(name,n,moist;rest=false)
    out=joinpath(root,name);mkpath(out)
    L=80000.0;H=40000.0;g=25.0;cp=13000.0;p0=5e5;T0=300.0
    constants=ThermodynamicConstants(Float64;gravitational_acceleration=g,dry_air_molar_mass=0.0023,dry_air_heat_capacity=cp)
    Rd=constants.molar_gas_constant/constants.dry_air.molar_mass;Rv=constants.molar_gas_constant/constants.vapor.molar_mass
    qbg=moist ? 1e-6 : 0.0
    grid=RectilinearGrid(GPU(),Float64;size=(n,n,n),x=(0,L),y=(-L/2,L/2),z=(0,H),halo=(4,4,4),topology=(Periodic,Periodic,Bounded))
    theta0=T0/(p0/1e5)^(Rd/cp)
    ref=ReferenceState(grid,constants;base_pressure=p0,potential_temperature=theta0,vapor_mass_fraction=qbg)
    micro=moist ? SaturationAdjustment(Float64;solver=SecantSolver(Float64;abstol=1e-7,maxiter=80)) : nothing
    f=2*1.7585e-4*sind(55.0)
    bounded=WENO(order=5,bounds=(0,1))
    scalars=moist ? (ρθ=WENO(order=5),ρqᵉ=bounded,dye=bounded) : (ρθ=WENO(order=5),ρqᵛ=bounded,dye=bounded)
    model=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref),thermodynamic_constants=constants,
        microphysics=micro,coriolis=FPlane(f=f),tracers=(:dye,),momentum_advection=WENO(order=5),scalar_advection=scalars,
        closure=ScalarDiffusivity(ν=50.0,κ=50.0))
    rho=vec(Array(interior(ref.density)));Tr=vec(Array(interior(ref.temperature)));p=vec(Array(interior(ref.pressure)))
    x=((1:n).-0.5).*L/n;y=((1:n).-0.5).*L/n.-L/2;z=((1:n).-0.5).*H/n
    # The two plume cases have matched initial density buoyancy, not simply identical T.
    bubble=[exp(-((xx-L/2)^2+yy^2+(zz-8000)^2)/(2*5000^2)) for xx in x,yy in y,zz in z]
    moisture=qbg .+ (moist && !rest ? 0.02 : 0.0).*bubble
    dT=rest ? 0.0 : 2.0
    Rref=(1-qbg)*Rd+qbg*Rv
    T=[Tr[k]*Rref/((1-moisture[i,j,k])*Rd+moisture[i,j,k]*Rv)*(1+dT*bubble[i,j,k]/Tr[k]) for i=1:n,j=1:n,k=1:n]
    set!(model;T,qᵗ=moisture,dye=reshape(rho,1,1,n).*bubble)
    sim=Simulation(model;Δt=0.25,stop_time=3600.0)
    conjure_time_step_wizard!(sim;cfl=0.2,max_Δt=3.0)
    Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
    heights=[4000.0,12000.0,24000.0];ks=[argmin(abs.(z.-h)) for h in heights]
    metadata=Dict("name"=>name,"nx"=>n,"ny"=>n,"nz"=>n,"lx"=>L,"ly"=>L,"height"=>H,"f"=>f,
        "meridional_boundary"=>"periodic","stop_seconds"=>3600.0,"output_seconds"=>120.0,"status"=>"RUNNING",
        "video_heights_requested_m"=>heights,"video_heights_actual_m"=>z[ks],"microphysics"=>moist ? "mixed-phase saturation adjustment" : "dry",
        "initialization"=>"Matched density buoyancy: 2 K dry-equivalent Gaussian at z=8 km; moist peak water mass fraction 0.020001",
        "scalar_transport"=>"WENO5; bounded moisture and dye (0,1); SSPRK3; summed advective CFL target 0.2",
        "rest_test"=>rest,"base_pressure_Pa"=>p0,"reference_bottom_T_K"=>T0,"gravity"=>g,"dry_molar_mass"=>0.0023,"dry_heat_capacity"=>cp,
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))),"breeze_commit"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b")
    open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,metadata);end
    history=NamedTuple[]
    function output(sim)
        u,v,w=map(fld->Array(interior(fld)),values(sim.model.velocities))
        uc=(u+circshift(u,(-1,0,0)))/2;vc=(v+circshift(v,(0,-1,0)))/2
        vort=(circshift(vc,(-1,0,0))-circshift(vc,(1,0,0)))/(2L/n)-(circshift(uc,(0,-1,0))-circshift(uc,(0,1,0)))/(2L/n)
        C=Array(interior(model.tracers.dye))./reshape(rho,1,1,n)
        qt=Array(interior(model.moisture_density))./reshape(rho,1,1,n)
        qv=moist ? Array(interior(model.microphysical_fields.qᵛ)) : qt
        ql=moist ? Array(interior(model.microphysical_fields.qˡ)) : zeros(size(qt))
        qi=moist ? Array(interior(model.microphysical_fields.qⁱ)) : zeros(size(qt))
        TT=Array(interior(model.temperature))
        @test all(isfinite,TT) && all(isfinite,w) && all(isfinite,qt)
        @test minimum(qt)>-1e-9
        @test minimum(C)>-1e-9 && maximum(C)<1+1e-9
        @test maximum(abs.(qv+ql+qi-qt))<1e-8
        water=(H/n)*sum(rho.*vec(mean(qt,dims=(1,2))))
        dye=(H/n)*sum(rho.*vec(mean(C,dims=(1,2))))
        wrms=sqrt(mean(abs2,w));cloud=maximum(ql+qi)
        push!(history,(t=sim.model.clock.time,water,dye,wrms,cloud,maxw=maximum(abs,w)))
        isapprox(water,history[1].water;rtol=1e-8,atol=1e-9) || error("Water inventory drift")
        isapprox(dye,history[1].dye;rtol=1e-8) || error("Passive tracer inventory drift")
        index=round(Int,sim.model.clock.time/120)
        jldsave(joinpath(out,@sprintf("frame_%05d.jld2",index));t=sim.model.clock.time,x,y,z,video_heights=z[ks],
            vorticity_slices=Float32.(vort[:,:,ks]),tracer_slices=Float32.(C[:,:,ks]),vapor_slices=Float32.(qv[:,:,ks]),
            condensate_slices=Float32.((ql+qi)[:,:,ks]),temperature_profile=vec(mean(TT,dims=(1,2))),
            temperature_xz=TT[:,n÷2,:],w_xz=w[:,n÷2,:],water,wrms,cloud)
        @printf("%s frame=%d t=%.1f maxw=%.4g cloud=%.4g\n",name,index,sim.model.clock.time,maximum(abs,w),cloud);flush(stdout)
    end
    add_callback!(sim,output,TimeInterval(120.0))
    run!(sim)
    rest && (@test maximum(h.wrms for h in history)<1e-4)
    moist && !rest && (@test maximum(h.cloud for h in history)>1e-6)
    metadata["status"]="COMPLETE";metadata["simulated_seconds"]=model.clock.time
    metadata["initial_max_condensate"]=history[1].cloud;metadata["peak_condensate"]=maximum(h.cloud for h in history)
    metadata["peak_w_rms"]=maximum(h.wrms for h in history)
    open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,metadata);end
    open(joinpath(out,"diagnostics.csv"),"w") do io
        println(io,join(string.(keys(history[1])),','));for row in history;println(io,join(values(row),','));end
    end
    return metadata
end
records=Dict{String,Any}()
for (name,n,moist,rest) in (("moist_rest",32,true,true),("dry_plume",64,false,false),("moist_plume",64,true,false),("moist_plume_fine",96,true,false))
    records[name]=plume(name,n,moist;rest)
    GC.gc();CUDA.reclaim()
end
records["status"]="PASS"
records["scope"]="Jovian-carrier implementation validation; no precipitation, ammonia chemistry, or calibrated Jovian moisture sounding. Resolution sensitivity must be assessed from the saved plume fields."
open(joinpath(root,"PASS.toml"),"w") do io;TOML.print(io,records);end

using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans.Units
using Oceananigans: interior
include("initial_turbulence.jl")

@inline function upper_sink_weight(z, bottom, top)
    s = (z-bottom)/(top-bottom)
    return 0 < s < 1 ? sinpi(s)^2 : zero(s)
end

@inline function tapered_diffusivity(z, p)
    s = clamp((z-(p.cap_height-p.width))/p.width, 0, 1)
    smooth = s*s*(3-2s)
    return p.lower + (p.cap-p.lower)*smooth
end

function main(args=ARGS)
length(args) >= 2 || error("usage: julia --project scripts/run.jl CONFIG OUTPUT_DIR [CPU|GPU]")
cfg = TOML.parsefile(args[1]); out = abspath(args[2])
isdir(out) && !isempty(readdir(out)) && error("Output exists and is nonempty: $out")
mkpath(out)
archname = length(args) >= 3 ? args[3] : "CPU"
archname in ("CPU", "GPU") || error("architecture must be CPU or GPU")
if archname == "GPU"
    CUDA.functional() || error("CUDA is not functional; refusing CPU fallback")
end
arch = archname == "GPU" ? GPU() : CPU()
FT = Float64 # establish a reference before evaluating Float32
Nx, Ny, Nz = cfg["nx"], cfg["ny"], cfg["nz"]
Lx, Ly, H = cfg["lx"], cfg["ly"], cfg["height"]
g, cp, theta0 = cfg["gravity"], cfg["heat_capacity"], cfg["theta0"]
zi, N2, Nc2 = cfg["cap_height"], cfg["cap_n2"], cfg["convective_n2"]
theta_ref(z) = theta0 * exp((Nc2 * min(z, zi) + N2 * max(z-zi, 0)) / g)
grid = RectilinearGrid(arch, FT; size=(Nx,Ny,Nz), x=(0,Lx), y=(0,Ly), z=(0,H), halo=(3,3,3), topology=(Periodic,Periodic,Bounded))
constants = ThermodynamicConstants(FT; gravitational_acceleration=g, dry_air_molar_mass=cfg["molar_mass"], dry_air_heat_capacity=cp)
ref = ReferenceState(grid, constants; base_pressure=cfg["base_pressure"], potential_temperature=theta_ref)
dynamics = AnelasticDynamics(ref)
f = 2cfg["rotation_rate"] * sind(cfg["latitude"])
z=collect(((1:Nz).-0.5).*H/Nz); dz=H/Nz
rho=vec(Array(interior(ref.density))); pressure=vec(Array(interior(ref.pressure)))
cooling_model=get(cfg,"cooling_model","relaxation")
cooling_model in ("relaxation","fixed_power") || error("Unknown cooling_model")
cooling_flux=get(cfg,"upper_cooling_flux",0.0)
cooling_flux >= 0 || error("upper_cooling_flux must be a nonnegative loss")
sink_weights=upper_sink_weight.(z,zi,H)
sink_normalization=sum(sink_weights)*dz
sink_normalization > 0 || error("Upper cooling layer must contain grid cells")
upper_energy_source=-cooling_flux.*sink_weights./sink_normalization
mixing_model=get(cfg,"mixing_model","uniform")
mixing_model in ("uniform","tapered") || error("Unknown mixing_model")
mixing=(; lower=cfg["thermal_diffusivity"], cap=get(cfg,"cap_diffusivity",cfg["thermal_diffusivity"]),
          cap_height=zi, width=get(cfg,"diffusivity_taper_width",20000.0))
0 < mixing.width <= zi || error("Invalid diffusivity taper width")
min(mixing.lower,mixing.cap) >= 0 || error("Diffusivity must be nonnegative")
kappa_profile=mixing_model=="tapered" ? tapered_diffusivity.(z,Ref(mixing)) : fill(mixing.lower,Nz)
kappa=mixing_model=="tapered" ? ((x,y,z,t)->tapered_diffusivity(z,mixing)) : mixing.lower
# theta is referenced to 1 bar: convert physical bottom enthalpy flux to rho*theta flux.
R = 8.314462618 / cfg["molar_mass"]
exner_bottom = (cfg["base_pressure"]/1e5)^(R/cp)
flux = cfg["heat_flux"] / (cp * exner_bottom)
bcs = (; ρθ=FieldBoundaryConditions(bottom=FluxBoundaryCondition(flux)))
p = (; theta0, zi, N2, Nc2, g, H, tau=cfg["cooling_time"], sponge=cfg["sponge_time"], U=cfg["jet_speed"], Ly, f)
@inline function cooling(x,y,z,t,theta,p)
    target = p.theta0 * exp((p.Nc2*min(z,p.zi)+p.N2*max(z-p.zi,0))/p.g)
    mask = clamp((z-p.zi)/(p.H-p.zi),0,1)
    return -mask*(theta-target)/p.tau
end
@inline damping(x,y,z,t,w,p) = -clamp((z-0.9p.H)/(0.1p.H),0,1)^2*w/p.sponge
# Optional prescribed geostrophic pressure-gradient surrogate, explicit in metadata.
@inline balance_v(x,y,z,t,p) = p.f*p.U*sin(2pi*y/p.Ly)
@inline fixed_sink(x,y,z,t,p) = -p.loss*upper_sink_weight(z,p.bottom,p.top)/p.normalization
thermal_forcing = if cooling_model=="fixed_power"
    # Volumetric enthalpy loss [W/m³]. Breeze converts this into a rho-theta tendency.
    (; ρE=Forcing(fixed_sink; parameters=(; loss=cooling_flux,bottom=zi,top=H,normalization=sink_normalization)))
else
    (; θ=Forcing(cooling; field_dependencies=:θ, parameters=p))
end
forcing = merge(thermal_forcing,(; w=Forcing(damping; field_dependencies=:w, parameters=p),
                                  v=Forcing(balance_v; parameters=p)))
model = AtmosphereModel(grid; dynamics, thermodynamic_constants=constants, coriolis=FPlane(f=f),
    advection=WENO(order=5), closure=ScalarDiffusivity(ν=cfg["viscosity"], κ=kappa),
    boundary_conditions=bcs, forcing)
rng=MersenneTwister(cfg["seed"])
noise=rand(rng,FT,Nx,Ny,Nz).-0.5
noise .-= mean(noise,dims=(1,2))
theta = [theta_ref((k-0.5)*H/Nz) + cfg["perturbation"]*noise[i,j,k]*exp(-((k-0.5)*H/Nz/(0.35H))^2) for i=1:Nx,j=1:Ny,k=1:Nz]
U = cfg["jet_speed"]
initial_velocity_info = Dict{String,Any}()
if get(cfg,"initial_velocity","rest") == "turbulence"
    U == 0 || error("Do not combine a prescribed jet with the turbulence initialization")
    u0,v0,w0,initial_velocity_info = initial_turbulence(cfg,rho,f)
    set!(model; θ=theta,u=u0,v=v0,w=w0)
else
    set!(model; θ=theta, u=(x,y,z)->U*sin(2pi*y/Ly))
end
thetafield=Field(liquid_ice_potential_temperature(model))
sim=Simulation(model; Δt=cfg["initial_dt"], stop_time=cfg["stop_seconds"])
conjure_time_step_wizard!(sim; cfl=0.4, max_Δt=cfg["max_dt"])
Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
metadata=merge(cfg,Dict("architecture"=>archname,"float_type"=>string(FT),"f"=>f,"julia_version"=>string(VERSION),"breeze_commit"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b","status"=>"RUNNING","physical_model"=>"dry anelastic, f-plane, periodic sides, specified bottom heat flux, free-slip impermeable boundaries, stable cap; no clouds"))
metadata["coriolis_model"]="f-plane"; metadata["beta"]=0.0
metadata["rotation_period_seconds"]=2pi/cfg["rotation_rate"]
metadata["script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
metadata["configuration_sha256"]=bytes2hex(sha256(read(args[1])))
metadata["cooling_model"]=cooling_model; metadata["mixing_model"]=mixing_model
metadata["discrete_fixed_upper_power_W_m2"]=cooling_model=="fixed_power" ? sum(upper_energy_source)*dz : 0.0
metadata["diagnostics_version"]=3
merge!(metadata,initial_velocity_info)
metadata["initial_turbulence_script_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"initial_turbulence.jl"))))
metadata["diagnostic_height_m"]=get(cfg,"diagnostic_height",zi)
open(joinpath(out,"metadata.toml"),"w") do io; TOML.print(io,metadata); end
theta_target=theta_ref.(z)
exner=(pressure./1e5).^(R/cp)
coolmask=clamp.((z.-zi)./(H-zi),0,1)
zf=(0:Nz).*dz
rhof=vcat(rho[1],(rho[1:end-1].+rho[2:end])./2,rho[end])
spongemask=clamp.((zf.-0.9H)./(0.1H),0,1).^2
diagnostic_seconds=min(get(cfg,"diagnostic_seconds",cfg["output_seconds"]),cfg["output_seconds"])
isapprox(cfg["output_seconds"]/diagnostic_seconds,round(cfg["output_seconds"]/diagnostic_seconds)) || error("output_seconds must be a multiple of diagnostic_seconds")
diagnostic_file=joinpath(out,"diagnostics.csv")
open(diagnostic_file,"w") do io
    println(io,"time_seconds,kinetic_energy_J_m2,w_rms_center_m_s,horizontal_rms_m_s,max_abs_w_m_s,theta_column_kg_K_m2,bottom_theta_flux_kg_K_m2_s,relaxation_theta_rate_kg_K_m2_s,relaxation_enthalpy_W_m2,sponge_power_W_m2,buoyancy_power_center_W_m2,vorticity_rms_over_f,theta_variance_K2,upper_theta_rate_kg_K_m2_s,upper_enthalpy_W_m2,external_enthalpy_net_W_m2,cap_contrast_ratio,cap_mean_profile_departure_rms_K")
end
started=time(); count=Ref(0)
timings=Tuple{Float64,Float64}[]
function save_frame(sim)
    m=sim.model; u,v,w=m.velocities
    compute!(thetafield)
    k=clamp(round(Int,get(cfg,"diagnostic_height",zi)/H*Nz),1,Nz)
    uc=Array(interior(u)); vc=Array(interior(v)); wc=Array(interior(w))
    th=Array(interior(thetafield))
    all(isfinite,th) && all(isfinite,wc) && all(isfinite,uc) && all(isfinite,vc) || error("non-finite state")
    # Face-weighted kinetic energy. Vertical boundary velocities are zero.
    ke=dz/2*(sum(rho .* vec(mean(uc.^2 .+ vc.^2,dims=(1,2)))) + sum(rhof .* vec(mean(wc.^2,dims=(1,2)))))
    sponge_power=-dz*sum(rhof .* spongemask .* vec(mean(wc.^2,dims=(1,2))))/cfg["sponge_time"]
    wcenter=(wc[:,:,1:end-1].+wc[:,:,2:end])./2
    theta_profile=vec(mean(th,dims=(1,2)))
    horizontal_ke_profile=rho.*vec(mean(uc.^2 .+vc.^2,dims=(1,2)))./2
    vertical_ke_profile=rho.*vec(mean(wcenter.^2,dims=(1,2)))./2
    theta_anomaly=th.-reshape(theta_profile,1,1,Nz)
    theta_column=dz*sum(rho .* theta_profile)
    # The legacy relaxation_* output names are retained as aliases for the total upper source.
    relaxation=cooling_model=="fixed_power" ? upper_energy_source./(rho.*cp.*exner) : -coolmask.*(theta_profile.-theta_target)./cfg["cooling_time"]
    relaxation_theta=dz*sum(rho.*relaxation)
    relaxation_enthalpy=dz*sum(rho.*cp.*exner.*relaxation)
    # Resolved transport and buoyancy work, interpolated to cell centers.
    # These are diagnostics, not the exact discrete WENO flux or a closed KE budget.
    heat_flux_profile=rho.*cp.*exner.*vec(mean(wcenter.*theta_anomaly,dims=(1,2)))
    buoyancy_power=dz*sum(rho.*g.*vec(mean(wcenter.*(1 .- reshape(theta_target,1,1,Nz)./th),dims=(1,2))))
    theta_gradient=vcat((theta_profile[2]-theta_profile[1])/dz,(theta_profile[3:end].-theta_profile[1:end-2])./(2dz),(theta_profile[end]-theta_profile[end-1])/dz)
    diffusive_heat_flux_profile=-rho.*cp.*exner.*kappa_profile.*theta_gradient
    mean_N2_profile=g.*theta_gradient./theta_profile
    cap_first=findfirst(>=(zi),z)
    cap_contrast_ratio=(theta_profile[end]-theta_profile[cap_first])/(theta_target[end]-theta_target[cap_first])
    cap_departure=sqrt(mean(abs2,theta_profile[cap_first:end].-theta_target[cap_first:end]))
    uc=(uc+circshift(uc,(-1,0,0)))/2; vc=(vc+circshift(vc,(0,-1,0)))/2
    zeta=(circshift(vc,(-1,0,0))-circshift(vc,(1,0,0)))/(2Lx/Nx) - (circshift(uc,(0,-1,0))-circshift(uc,(0,1,0)))/(2Ly/Ny)
    values=(m.clock.time,ke,sqrt(mean(abs2,wcenter)),sqrt(mean(uc.^2 .+ vc.^2)),maximum(abs,wc),theta_column,flux,relaxation_theta,relaxation_enthalpy,sponge_power,buoyancy_power,sqrt(mean(abs2,zeta))/abs(f),mean(abs2,theta_anomaly),relaxation_theta,relaxation_enthalpy,cfg["heat_flux"]+relaxation_enthalpy,cap_contrast_ratio,cap_departure)
    all(isfinite,values) || error("non-finite diagnostic")
    open(diagnostic_file,"a") do io; println(io,join(values,",")); end
    save_output=isapprox(m.clock.time/cfg["output_seconds"],round(m.clock.time/cfg["output_seconds"]);atol=1e-8) || m.clock.time >= cfg["stop_seconds"]
    save_output || return nothing
    file=joinpath(out,@sprintf("frame_%05d.jld2",count[]))
    jldsave(file; t=m.clock.time, u=uc[:,:,k],v=vc[:,:,k],w=wc[:,:,k],zeta=zeta[:,:,k],
       theta=th[:,:,k],theta_xz=th[:,Ny÷2,:],w_xz=wc[:,Ny÷2,:],
       theta_profile,horizontal_ke_profile,vertical_ke_profile,w_rms=sqrt(mean(abs2,wc)),heat_flux_profile,diffusive_heat_flux_profile,kappa_profile,upper_energy_source,mean_N2_profile,
       kinetic_energy=ke,relaxation_theta,relaxation_enthalpy,sponge_power,buoyancy_power,
       x=((1:Nx).-0.5).*Lx/Nx,z,pressure,density=rho)
    @printf("frame=%d t=%.1fs iter=%d dt=%.3fs maxw=%.4g wall=%.1fs\n",count[],m.clock.time,m.clock.iteration,sim.Δt,maximum(abs,wc),time()-started)
    push!(timings,(Float64(m.clock.time),time()-started))
    flush(stdout);count[]+=1
end
add_callback!(sim,save_frame,TimeInterval(diagnostic_seconds))
run!(sim)
metadata["status"]="COMPLETE";metadata["wall_seconds"]=time()-started
metadata["iterations"]=model.clock.iteration;metadata["simulated_seconds"]=model.clock.time
if length(timings) >= 3
    metadata["post_warmup_wall_seconds"]=last(timings)[2]-timings[2][2]
    metadata["post_warmup_simulated_seconds"]=last(timings)[1]-timings[2][1]
end
open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,metadata);end

end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

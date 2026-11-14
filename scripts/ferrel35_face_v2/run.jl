using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf, SHA
using Oceananigans.Units
using Oceananigans: interior
include("initial_flow.jl")
include("zonal_output.jl")
include("native_budget.jl")
include("face_budget.jl")
include("centered_interpretation.jl")
include("thermal_budget_controller.jl")

# The pinned Simulation checkpoint omits sim.Δt. Carry it in callback state.
mutable struct SavedTimeStep
    value::Float64
end
(s::SavedTimeStep)(sim) = (s.value=sim.Δt; nothing)
Oceananigans.prognostic_state(s::SavedTimeStep) = (; value=s.value)
Oceananigans.restore_prognostic_state!(s::SavedTimeStep, state::NamedTuple) = (s.value=state.value; s)

@inline function upper_sink_weight(z, bottom, top)
    s = (z-bottom)/(top-bottom)
    return 0 < s < 1 ? sinpi(s)^2 : zero(s)
end

@inline function tapered_diffusivity(z, p)
    s = clamp((z-(p.cap_height-p.width))/p.width, 0, 1)
    smooth = s*s*(3-2s)
    return p.lower + (p.cap-p.lower)*smooth
end

function run_eddy(args=ARGS)
length(args) >= 2 || error("usage: julia --project scripts/eddy1200/run.jl CONFIG OUTPUT_DIR [CPU|GPU] [CHECKPOINT]")
cfg = TOML.parsefile(args[1]); out = abspath(args[2])
pickup = length(args) >= 4 ? abspath(args[4]) : false
pickup === false || isfile(pickup) || error("Missing restart checkpoint")
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
get(cfg,"meridional_boundary","periodic") in ("periodic","channel") || error("Unknown meridional boundary")
get(cfg,"rotation_model","f-plane") in ("f-plane","beta-plane") || error("Unknown rotation model")
channel = get(cfg,"meridional_boundary","periodic") == "channel"
beta = get(cfg,"beta",0.0)
get(cfg,"rotation_model","f-plane") == "f-plane" && beta != 0 && error("F-plane requires beta=0")
!channel && beta != 0 && error("Linear beta requires a bounded meridional channel")
Ytopology = channel ? Bounded : Periodic
ybounds = channel ? (-Ly/2,Ly/2) : (0,Ly)
grid = RectilinearGrid(arch, FT; size=(Nx,Ny,Nz), x=(0,Lx), y=ybounds, z=(0,H), halo=(4,4,4), topology=(Periodic,Ytopology,Bounded))
constants = ThermodynamicConstants(FT; gravitational_acceleration=g, dry_air_molar_mass=cfg["molar_mass"], dry_air_heat_capacity=cp)
ref = ReferenceState(grid, constants; base_pressure=cfg["base_pressure"], potential_temperature=theta_ref)
dynamics = AnelasticDynamics(ref)
f = 2cfg["rotation_rate"] * sind(cfg["latitude"])
z=collect(((1:Nz).-0.5).*H/Nz); dz=H/Nz
rho=vec(Array(interior(ref.density))); pressure=vec(Array(interior(ref.pressure)))
cooling_model=get(cfg,"cooling_model","stable_cap")
cooling_model=="stable_cap" || error("Shear experiment requires stable_cap")
cooling_model in ("relaxation","fixed_power","stable_cap") || error("Unknown cooling_model")
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
cap_width=get(cfg,"cap_relaxation_transition",10000.0)
0 < cap_width < H-zi || error("Invalid cap relaxation transition")
rate=archname=="GPU" ? CUDA.zeros(FT,1) : zeros(FT,1)
p = (; rate, theta0, zi, N2, Nc2, g, H, cap_width, stable_cap=cooling_model=="stable_cap", tau=cfg["cooling_time"], sponge=cfg["sponge_time"], U=cfg["jet_speed"], Ly, f)
@inline function cooling(x,y,z,t,theta,p)
    target = p.theta0 * exp((p.Nc2*min(z,p.zi)+p.N2*max(z-p.zi,0))/p.g)
    s = clamp((z-p.zi)/(p.stable_cap ? p.cap_width : p.H-p.zi),0,1)
    mask = p.stable_cap ? s*s*(3-2s) : s
    return -mask*(theta-target)/p.tau - p.rate[1]*theta
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
model = AtmosphereModel(grid; dynamics, thermodynamic_constants=constants, coriolis=(get(cfg,"rotation_model","f-plane")=="beta-plane" ? BetaPlane(f₀=f, β=beta) : FPlane(f=f)),
    tracers=(:dye,), momentum_advection=WENO(order=5),
    scalar_advection=(ρθ=WENO(order=5), ρqᵛ=WENO(order=5,bounds=(0,1)), dye=WENO(order=5,bounds=(0,1))), closure=ScalarDiffusivity(ν=cfg["viscosity"], κ=(ρθ=kappa, ρqᵛ=kappa, dye=cfg["tracer_diffusivity"])),
    boundary_conditions=bcs, forcing)
rng=MersenneTwister(cfg["seed"])
noise=rand(rng,FT,Nx,Ny,Nz).-0.5
noise .-= mean(noise,dims=(1,2))
theta = [theta_ref((k-0.5)*H/Nz) + cfg["perturbation"]*noise[i,j,k]*exp(-((k-0.5)*H/Nz/(0.35H))^2) for i=1:Nx,j=1:Ny,k=1:Nz]
U = cfg["jet_speed"]
U==0 || error("Legacy imposed jet/body force must be disabled")
S=cfg["shear_rossby"]*f
Us=S*Ly/pi
background_u(y)=-Us*sinpi(y/Ly)
# ∂y(p′/rho_ref) = −f(y) U(y); independent of z for this barotropic flow.
k_s=pi/Ly
background_pressure(y)=Us*(-f*cos(k_s*y)/k_s + beta*(-y*cos(k_s*y)/k_s + sin(k_s*y)/k_s^2))
initial_velocity_info = Dict{String,Any}()
if get(cfg,"initial_velocity","rest") == "turbulence"
    U == 0 || error("Do not combine a prescribed jet with the turbulence initialization")
    channel || error("Large-eddy initialization requires channel walls")
    u0,v0,w0,p0,initial_velocity_info = initial_large_eddies(cfg,rho,f)
    yc=((1:Ny).-0.5).*Ly/Ny .+ first(ybounds)
    u0 .+= reshape(background_u.(yc),1,Ny,1)
    set!(model; θ=theta,u=u0,v=v0,w=w0)
else
    set!(model; θ=theta, u=(x,y,z)->background_u(y))
end
if @isdefined(p0)
    p0 .+= reshape(background_pressure.(yc),1,Ny,1)
    set!(model.dynamics.pressure_anomaly,p0)
    p0=nothing
else
    set!(model.dynamics.pressure_anomaly,(x,y,z)->background_pressure(y))
end
xc=((1:Nx).-0.5).*Lx/Nx;yc=((1:Ny).-0.5).*Ly/Ny .+ first(ybounds)
function tracer_patch(x,y)
    get(cfg,"tracer_uniform",false) && return 0.5
    sigma=0.08min(Lx,Ly)
    blob(x0,y0)=exp(-0.5*((Lx/pi*sinpi((x-x0)/Lx)/sigma)^2+((y-y0)/sigma)^2))
    return 0.8blob(0.25Lx,-0.15Ly)+0.6blob(0.55Lx,0.1Ly)+0.5blob(0.8Lx,-0.25Ly)
end
concentration=[tracer_patch(x,y) for x in xc,y in yc]
dye0=reshape(concentration,Nx,Ny,1).*reshape(rho,1,1,Nz)
# AtmosphereModel user tracers are stored as tracer mass density between steps.
set!(model.tracers.dye,dye0)
dye0=nothing;concentration=nothing
noise=nothing; theta=nothing
if @isdefined(u0); u0=nothing; v0=nothing; w0=nothing; end
GC.gc()
thetafield=Field(liquid_ice_potential_temperature(model))
sim=Simulation(model; Δt=cfg["initial_dt"], stop_time=cfg["stop_seconds"])
conjure_time_step_wizard!(sim; cfl=0.2, max_Δt=cfg["max_dt"])
Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
saved_dt=SavedTimeStep(sim.Δt)
add_callback!(sim,saved_dt,IterationInterval(1);name=:saved_time_step)
metadata=merge(cfg,Dict("architecture"=>archname,"float_type"=>string(FT),"f"=>f,"julia_version"=>string(VERSION),"breeze_commit"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b","status"=>"RUNNING","physical_model"=>"dry anelastic, f-plane, periodic sides, specified bottom heat flux, free-slip impermeable boundaries, stable cap; no clouds"))
metadata["coriolis_model"]=get(cfg,"rotation_model","f-plane"); metadata["beta"]=beta
metadata["physical_model"]="dry anelastic, configured rotation and meridional boundaries, periodic zonal boundary, bottom heat flux, free-slip impermeable walls; no clouds"
metadata["meridional_boundary"]=channel ? "channel" : "periodic"
metadata["f_south"]=f-beta*Ly/2; metadata["f_north"]=f+beta*Ly/2
metadata["initial_channel_script_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"initial_flow.jl"))))
metadata["rotation_period_seconds"]=2pi/cfg["rotation_rate"]
metadata["script_sha256"]=bytes2hex(sha256(read(@__FILE__)))
metadata["configuration_sha256"]=bytes2hex(sha256(read(args[1])))
metadata["cooling_model"]=cooling_model; metadata["mixing_model"]=mixing_model
metadata["discrete_fixed_upper_power_W_m2"]=cooling_model=="fixed_power" ? sum(upper_energy_source)*dz : 0.0
metadata["diagnostics_version"]=7
metadata["scalar_transport"]="WENO5; bounds (0,1) on vapor and dye; SSPRK3; summed advective CFL target 0.2"
metadata["background_max_speed_m_s"]=abs(Us)
metadata["background_peak_zeta_over_f"]=cfg["shear_rossby"]
metadata["shear_initialization"]="U=-Us*sin(pi*y/Ly); balanced kinematic pressure; no maintained shear forcing"
metadata["tracer_definition"]="dimensionless passive mass fraction; initialized as three smooth patches; stored prognostically as rho*C"
metadata["controller_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"thermal_budget_controller.jl"))))
metadata["gpu_name"]=archname=="GPU" ? CUDA.name(CUDA.device()) : "none"
metadata["pickup_checkpoint"]=pickup === false ? "" : pickup
metadata["cap_treatment"]="Local cap relaxation plus globally uniform fractional theta correction balancing measured net external enthalpy input"
metadata["video_heights_requested_m"]=cfg["video_heights"]
video_indices=[argmin(abs.(z.-h)) for h in cfg["video_heights"]]
length(unique(video_indices)) == 3 || error("Require three distinct video depths")
metadata["video_heights_actual_m"]=z[video_indices]
metadata["vorticity_video_storage"]="Float32 slices; dynamics and bulk diagnostics Float64"
merge!(metadata,initial_velocity_info)
metadata["initial_turbulence_script_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"initial_flow.jl"))))
metadata["diagnostic_height_m"]=get(cfg,"diagnostic_height",zi)
open(joinpath(out,"metadata.toml"),"w") do io; TOML.print(io,metadata); end
theta_target=theta_ref.(z)
exner=(pressure./1e5).^(R/cp)
coolmask=clamp.((z.-zi)./(cooling_model=="stable_cap" ? cap_width : H-zi),0,1)
cooling_model=="stable_cap" && (coolmask=coolmask.^2 .* (3 .- 2coolmask))
cap_core=(z .>= get(cfg,"cap_core_bottom",90000.0)) .& (z .<= get(cfg,"cap_core_top",115000.0))
any(cap_core) || error("No grid cells in protected cap core")
theta_mean_field=Field(Average(thetafield,dims=(1,2)))
controller=ThermalBudgetController(rate,thetafield,theta_mean_field,rho,rho.*cp.*exner,coolmask,theta_target,
    dz,cfg["cooling_time"],cfg["heat_flux"],flux,get(cfg,"balance_heat",true),0.0,-1.0,0.0,0.0,0.0,0.0,0.0)
controller(sim)
add_callback!(sim,controller,TimeInterval(cfg["controller_seconds"]);name=:thermal_budget)
zf=(0:Nz).*dz
rhof=vcat(rho[1],(rho[1:end-1].+rho[2:end])./2,rho[end])
spongemask=clamp.((zf.-0.9H)./(0.1H),0,1).^2
diagnostic_seconds=min(get(cfg,"diagnostic_seconds",cfg["output_seconds"]),cfg["output_seconds"])
isapprox(cfg["output_seconds"]/diagnostic_seconds,round(cfg["output_seconds"]/diagnostic_seconds)) || error("output_seconds must be a multiple of diagnostic_seconds")
diagnostic_file=joinpath(out,"diagnostics.csv")
open(diagnostic_file,"w") do io
    println(io,"time_seconds,kinetic_energy_J_m2,w_rms_center_m_s,horizontal_rms_m_s,max_abs_w_m_s,theta_column_kg_K_m2,bottom_theta_flux_kg_K_m2_s,relaxation_theta_rate_kg_K_m2_s,relaxation_enthalpy_W_m2,sponge_power_W_m2,buoyancy_power_center_W_m2,vorticity_rms_over_f,theta_variance_K2,upper_theta_rate_kg_K_m2_s,upper_enthalpy_W_m2,external_enthalpy_net_W_m2,cap_contrast_ratio,cap_mean_profile_departure_rms_K,cap_core_min_N2_s2,cap_core_mean_N2_s2,tracer_column_kg_m2,tracer_min,tracer_max,cap_enthalpy_W_m2,budget_correction_W_m2,budget_alpha_s_inv,preupdate_net_W_m2,integrated_theta_source_kg_K_m2,integrated_external_enthalpy_J_m2,max_preupdate_net_W_m2")
end
native_budget=NativeBudget(model,out,cfg)
native_budget.face=FaceBudget(model,cfg,rate)
add_callback!(sim,native_budget,IterationInterval(1);name=:native_momentum_budget,callsite=TendencyCallsite())
started=time(); count=Ref(0)
# Keep checkpoints outside the frame directory; discard only superseded checkpoints.
checkpoint_dir=joinpath(out,"checkpoints")
sim.output_writers[:checkpointer]=Checkpointer(model; schedule=TimeInterval(cfg["checkpoint_seconds"]), dir=checkpoint_dir, cleanup=true, verbose=true)
timings=Tuple{Float64,Float64}[]
function save_frame(sim)
    m=sim.model; u,v,w=m.velocities
    compute!(thetafield)
    k=video_indices[1]
    uc=Array(interior(u)); vc=Array(interior(v)); wc=Array(interior(w))
    # Bounded meridional velocity has Ny+1 faces. Average squared faces for KE.
    vface = vc
    vsquared = channel ? (vc[:,1:end-1,:].^2 .+ vc[:,2:end,:].^2)./2 : vc.^2
    th=Array(interior(thetafield))
    dye_density=Array(interior(m.tracers.dye))
    dye=dye_density./reshape(rho,1,1,Nz)
    tracer_column=dz*sum(vec(mean(dye_density,dims=(1,2))))
    tracer_min=minimum(dye);tracer_max=maximum(dye)
    all(isfinite,dye) || error("Nonfinite passive tracer")
    all(isfinite,th) && all(isfinite,wc) && all(isfinite,uc) && all(isfinite,vc) || error("non-finite state")
    # Face-weighted kinetic energy. Vertical boundary velocities are zero.
    ke=dz/2*(sum(rho .* vec(mean(uc.^2 .+ vsquared,dims=(1,2)))) + sum(rhof .* vec(mean(wc.^2,dims=(1,2)))))
    sponge_power=-dz*sum(rhof .* spongemask .* vec(mean(wc.^2,dims=(1,2))))/cfg["sponge_time"]
    wcenter=(wc[:,:,1:end-1].+wc[:,:,2:end])./2
    theta_profile=vec(mean(th,dims=(1,2)))
    horizontal_ke_profile=rho.*vec(mean(uc.^2 .+vsquared,dims=(1,2)))./2
    vertical_ke_profile=rho.*vec(mean(wcenter.^2,dims=(1,2)))./2
    theta_anomaly=th.-reshape(theta_profile,1,1,Nz)
    theta_column=dz*sum(rho .* theta_profile)
    # The legacy relaxation_* output names are retained as aliases for the total upper source.
    relaxation=cooling_model=="fixed_power" ? upper_energy_source./(rho.*cp.*exner) : -coolmask.*(theta_profile.-theta_target)./cfg["cooling_time"]
    cap_enthalpy=dz*sum(rho.*cp.*exner.*relaxation)
    budget_correction=-controller.alpha*dz*sum(rho.*cp.*exner.*theta_profile)
    relaxation .-= controller.alpha.*theta_profile
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
    cap_min=minimum(mean_N2_profile[cap_core]); cap_mean=mean(mean_N2_profile[cap_core])
    if cooling_model=="stable_cap" && cap_min < cfg["cap_min_fraction"]*N2
        error("Protected cap core lost stability: min N2=$cap_min at t=$(m.clock.time)")
    end
    uc=(uc+circshift(uc,(-1,0,0)))/2; vc=channel ? (vc[:,1:end-1,:]+vc[:,2:end,:])/2 : (vc+circshift(vc,(0,-1,0)))/2
    du_dy = if channel
        # Even tangential extension through free-slip walls; no periodic seam.
        padded=cat(uc[:,1:1,:],uc,uc[:,end:end,:];dims=2)
        (padded[:,3:end,:].-padded[:,1:end-2,:])./(2Ly/Ny)
    else
        (circshift(uc,(0,-1,0))-circshift(uc,(0,1,0)))/(2Ly/Ny)
    end
    zeta=(circshift(vc,(-1,0,0))-circshift(vc,(1,0,0)))/(2Lx/Nx) - du_dy
    values=(m.clock.time,ke,sqrt(mean(abs2,wcenter)),sqrt(mean(uc.^2 .+ vc.^2)),maximum(abs,wc),theta_column,flux,relaxation_theta,relaxation_enthalpy,sponge_power,buoyancy_power,sqrt(mean(abs2,zeta))/abs(f),mean(abs2,theta_anomaly),relaxation_theta,relaxation_enthalpy,cfg["heat_flux"]+relaxation_enthalpy,cap_contrast_ratio,cap_departure,cap_min,cap_mean,tracer_column,tracer_min,tracer_max,cap_enthalpy,budget_correction,controller.alpha,controller.preupdate_net,controller.theta_integral,controller.net_enthalpy_integral,controller.max_preupdate_net)
    all(isfinite,values) || error("non-finite diagnostic")
    open(diagnostic_file,"a") do io; println(io,join(values,",")); end
    save_output=isapprox(m.clock.time/cfg["output_seconds"],round(m.clock.time/cfg["output_seconds"]);atol=1e-8) || m.clock.time >= cfg["stop_seconds"]
    save_output || return nothing
    count[]=round(Int,m.clock.time/cfg["output_seconds"])
    file=joinpath(out,@sprintf("frame_%05d.jld2",count[]))
    isfile(file) && error("Refusing to overwrite a saved frame")
    weights=reshape(rho./sum(rho),1,1,Nz)
    ubar=dropdims(sum(weights.*uc;dims=3);dims=3)
    vbar=dropdims(sum(weights.*vc;dims=3);dims=3)
    Rxx=dropdims(sum(weights.*uc.^2;dims=3);dims=3).-ubar.^2
    Rxy=dropdims(sum(weights.*uc.*vc;dims=3);dims=3).-ubar.*vbar
    Ryy=dropdims(sum(weights.*vc.^2;dims=3);dims=3).-vbar.^2
    wall_normal_max=channel ? maximum(abs,vface[:,[1,end],:]) : 0.0
    wall_normal_max < 1e-10 || error("Nonzero meridional wall velocity")
    zonal = zonal_output(m, uc, vc, wcenter, th, dye, Ly/Ny, dz)
    cumulative=cumulative_native_budget(native_budget,m)
    face_report=face_snapshot(native_budget.face,m)
    native_rhs=native_terms(native_budget)
    centered=centered_flux_interpretation(zonal.zonal_u,zonal.zonal_v,zonal.zonal_w,
        zonal.eddy_uv.*reshape(rho,1,:),zonal.eddy_uw.*reshape(rho,1,:),rho,Ly/Ny,dz,native_rhs.parts[:,:,1])
    jldsave(file; zonal..., cumulative..., face_report..., centered..., native_rhou_rhs=native_rhs.parts, ubar,vbar,Rxx,Rxy,Ryy,wall_normal_max,
       tracer_slices=Float32.(dye[:,:,video_indices]),tracer_column,tracer_min,tracer_max,cap_enthalpy,budget_correction,
       vorticity_slices=Float32.(zeta[:,:,video_indices]),video_heights=z[video_indices],cap_min,cap_mean,
       y=(collect(1:Ny).-0.5).*Ly/Ny .+ first(ybounds), t=m.clock.time, u=uc[:,:,k],v=vc[:,:,k],w=wc[:,:,k],zeta=zeta[:,:,k],
       theta=th[:,:,k],theta_xz=th[:,Ny÷2,:],w_xz=wc[:,Ny÷2,:],
       theta_profile,horizontal_ke_profile,vertical_ke_profile,w_rms=sqrt(mean(abs2,wc)),heat_flux_profile,diffusive_heat_flux_profile,kappa_profile,upper_energy_source,mean_N2_profile,
       kinetic_energy=ke,relaxation_theta,relaxation_enthalpy,sponge_power,buoyancy_power,
       x=((1:Nx).-0.5).*Lx/Nx,z,pressure,density=rho)
    @printf("frame=%d t=%.1fs iter=%d dt=%.3fs maxw=%.4g wall=%.1fs\n",count[],m.clock.time,m.clock.iteration,sim.Δt,maximum(abs,wc),time()-started)
    push!(timings,(Float64(m.clock.time),time()-started))
    flush(stdout);count[]+=1
end
add_callback!(sim,save_frame,TimeInterval(diagnostic_seconds))
# Segment output is separate, so a restored callback never overwrites earlier frames.
if pickup !== false
    set!(sim; checkpoint=pickup)
    sim.Δt=saved_dt.value
end
initialize_pressure_diagnostic!(model)
metadata["face_budget_script_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"face_budget.jl"))))
metadata["native_budget_scope"]="Every complete SSPRK3 step; cumulative native zonal momentum terms; sparse stage files; centered eddy decomposition remainder retained"
metadata["native_budget_script_sha256"]=bytes2hex(sha256(read(joinpath(@__DIR__,"native_budget.jl"))))
metadata["segment_start_seconds"]=model.clock.time
open(joinpath(out,"metadata.toml"),"w") do io; TOML.print(io,metadata); end
if pickup !== false
    # Simulation initialization refreshes the first RHS without invoking model
    # callbacks after restart. Seed this observer explicitly before that step.
    # update_state! refreshes auxiliary fields/RHS without advancing the state.
    Oceananigans.TimeSteppers.update_state!(model;compute_tendencies=true)
    model.clock.stage==1 || error("Restart is not at a full-step boundary")
    native_budget(model)
end
run!(sim; checkpoint_at_end=true)
flush_native_budget!(native_budget,model)
native_summary=cumulative_native_budget(native_budget,model)
jldsave(joinpath(out,"face_budget_final.jld2"); face_snapshot(native_budget.face,model)...)
native_budget.checks>0 || error("No complete native budget steps")
metadata["native_budget_checks"]=native_budget.checks
metadata["native_budget_max_absolute_residual"]=native_budget.max_absolute_residual
metadata["native_budget_cumulative_relative_residual"]=native_summary.native_integrated_relative_residual
metadata["status"]="COMPLETE";metadata["wall_seconds"]=time()-started
metadata["final_dt"]=sim.Δt
metadata["iterations"]=model.clock.iteration;metadata["simulated_seconds"]=model.clock.time
if length(timings) >= 3
    metadata["post_warmup_wall_seconds"]=last(timings)[2]-timings[2][2]
    metadata["post_warmup_simulated_seconds"]=last(timings)[1]-timings[2][1]
end
open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,metadata);end

end
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_eddy()
end

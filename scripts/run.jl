using Breeze, Oceananigans, CUDA, JLD2, TOML, Random, Statistics, Printf
using Oceananigans.Units
using Oceananigans: interior

function main()
length(ARGS) >= 2 || error("usage: julia --project scripts/run.jl CONFIG OUTPUT_DIR [CPU|GPU]")
cfg = TOML.parsefile(ARGS[1]); out = abspath(ARGS[2])
isdir(out) && !isempty(readdir(out)) && error("Output exists and is nonempty: $out")
mkpath(out)
archname = length(ARGS) >= 3 ? ARGS[3] : "CPU"
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
forcing = (; θ=Forcing(cooling; field_dependencies=:θ, parameters=p),
             w=Forcing(damping; field_dependencies=:w, parameters=p),
             v=Forcing(balance_v; parameters=p))
model = AtmosphereModel(grid; dynamics, thermodynamic_constants=constants, coriolis=FPlane(f=f),
    advection=WENO(order=5), closure=ScalarDiffusivity(ν=cfg["viscosity"], κ=cfg["thermal_diffusivity"]),
    boundary_conditions=bcs, forcing)
rng=MersenneTwister(cfg["seed"])
noise=rand(rng,FT,Nx,Ny,Nz).-0.5
noise .-= mean(noise,dims=(1,2))
theta = [theta_ref((k-0.5)*H/Nz) + cfg["perturbation"]*noise[i,j,k]*exp(-((k-0.5)*H/Nz/(0.35H))^2) for i=1:Nx,j=1:Ny,k=1:Nz]
U = cfg["jet_speed"]
set!(model; θ=theta, u=(x,y,z)->U*sin(2pi*y/Ly))
thetafield=Field(liquid_ice_potential_temperature(model))
sim=Simulation(model; Δt=cfg["initial_dt"], stop_time=cfg["stop_seconds"])
conjure_time_step_wizard!(sim; cfl=0.4, max_Δt=cfg["max_dt"])
Oceananigans.Diagnostics.erroring_NaNChecker!(sim)
metadata=merge(cfg,Dict("architecture"=>archname,"float_type"=>string(FT),"f"=>f,"julia_version"=>string(VERSION),"breeze_commit"=>"0479d9c8ce006adcb1e3314c324b84980ad4871b","status"=>"RUNNING","physical_model"=>"dry anelastic, f-plane, periodic sides, heated free-slip bottom, stable cap; no clouds"))
open(joinpath(out,"metadata.toml"),"w") do io; TOML.print(io,metadata); end
started=time(); count=Ref(0)
function save_frame(sim)
    m=sim.model; u,v,w=m.velocities
    compute!(thetafield)
    k=clamp(round(Int,zi/H*Nz),1,Nz)
    uc=Array(interior(u)); vc=Array(interior(v)); wc=Array(interior(w))
    th=Array(interior(thetafield))
    all(isfinite,th) && all(isfinite,wc) || error("non-finite state")
    uc=(uc+circshift(uc,(-1,0,0)))/2; vc=(vc+circshift(vc,(0,-1,0)))/2
    zeta=(circshift(vc,(-1,0,0))-circshift(vc,(1,0,0)))/(2Lx/Nx) - (circshift(uc,(0,-1,0))-circshift(uc,(0,1,0)))/(2Ly/Ny)
    file=joinpath(out,@sprintf("frame_%05d.jld2",count[]))
    jldsave(file; t=m.clock.time, u=uc[:,:,k],v=vc[:,:,k],w=wc[:,:,k],zeta=zeta[:,:,k],
       theta=th[:,:,k],theta_xz=th[:,Ny÷2,:],w_xz=wc[:,Ny÷2,:],
       theta_profile=vec(mean(th,dims=(1,2))),w_rms=sqrt(mean(abs2,wc)),
       x=((1:Nx).-0.5).*Lx/Nx, z=((1:Nz).-0.5).*H/Nz,
       pressure=vec(Array(interior(ref.pressure))), density=vec(Array(interior(ref.density))))
    @printf("frame=%d t=%.1fs iter=%d dt=%.3fs maxw=%.4g wall=%.1fs\n",count[],m.clock.time,m.clock.iteration,sim.Δt,maximum(abs,wc),time()-started)
    flush(stdout);count[]+=1
end
add_callback!(sim,save_frame,TimeInterval(cfg["output_seconds"]))
run!(sim)
metadata["status"]="COMPLETE";metadata["wall_seconds"]=time()-started
metadata["iterations"]=model.clock.iteration;metadata["simulated_seconds"]=model.clock.time
open(joinpath(out,"metadata.toml"),"w") do io;TOML.print(io,metadata);end

end
main()

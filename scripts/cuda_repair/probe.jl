using CUDA, Oceananigans, Breeze, TOML, SHA, Dates
CUDA.functional() || error("GPU runtime unavailable")
CUDA.versioninfo()
@assert v"12.2" <= CUDA.runtime_version() < v"12.3"
x=CUDA.fill(2.0,1024); x .= 3 .* x .+ 1
CUDA.synchronize()
@assert sum(x)==7168.0
# Exercise the GPU anelastic pressure solver and time stepping.
grid=RectilinearGrid(GPU();size=(16,16,16),extent=(1000,1000,1000),halo=(4,4,4))
c=ThermodynamicConstants(Float64; gravitational_acceleration=25.0,dry_air_molar_mass=0.0023,dry_air_heat_capacity=13000.0)
ref=ReferenceState(grid,c;base_pressure=2e6,potential_temperature=170.0)
m=AtmosphereModel(grid;dynamics=AnelasticDynamics(ref),thermodynamic_constants=c,coriolis=FPlane(f=2e-4))
set!(m;θ=(x,y,z)->170.0+0.01*sin(2pi*x/1000)*sin(2pi*y/1000)*sin(pi*z/1000))
s=Simulation(m;Δt=0.01,stop_iteration=3);run!(s)
CUDA.synchronize()
@assert m.clock.iteration==3
out=only(ARGS);mkpath(out)
open(joinpath(out,"PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","utc"=>string(now(UTC)),"runtime"=>string(CUDA.runtime_version()),"gpu"=>CUDA.name(CUDA.device()),"iterations"=>m.clock.iteration,"hashes"=>Dict(f=>bytes2hex(sha256(read(f))) for f in ["Project.toml","Manifest.toml","environments/cascade-cuda122/Project.toml","scripts/cuda_repair/environment.sh","scripts/cuda_repair/probe.jl"])))
end

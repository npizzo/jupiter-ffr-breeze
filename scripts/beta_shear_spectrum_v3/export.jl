include("../beta_shear_repair_v5/model.jl")
include("../beta_shear_repair_v5/initialization.jl")
include("../beta_shear_repair_v5/cap_initialization.jl")
function main()
    out=ARGS[1];mkpath(out)
    base=merge(TOML.parsefile(joinpath(@__DIR__,"../beta_shear_repair_v5/pilot.toml")),Dict("nx"=>8,"lx"=>.25,"initial_rms"=>0.,"perturbation"=>0.,"heat"=>0.))
    for (ny,nz) in [(128,80),(256,80)], (name,opts) in [
        ("reference",Dict("jet_speed"=>0.)),
        ("deep_anti",Dict("jet_speed"=>.25)),
        ("cap_anti",Dict("jet_speed"=>.015625,"jet_design"=>"cap")),
        ("cap_cyclone",Dict("jet_speed"=>-.015625,"jet_design"=>"cap")),
        ("traditional_anti",Dict("jet_speed"=>.25,"horizontal_rotation_factor"=>0.)),
        ("traditional_cyclone",Dict("jet_speed"=>-.25,"horizontal_rotation_factor"=>0.))]
        cfg=merge(base,opts,Dict("ny"=>ny,"nz"=>nz));m,p=build_model(cfg;arch=CPU());a=arrays(m)
        file=joinpath(out,"$(name)_$(ny)_$(nz).jld2")
        jldsave(file;rho=a.rho[1,:,:],T=a.T[1,:,:],u=a.uf[1,:,:],dy=p.dy,dz=p.dz,f=p.f,h=p.fh,beta=p.beta,ly=p.ly)
        println("EXPORTED ",file);flush(stdout);GC.gc()
    end
end
main()

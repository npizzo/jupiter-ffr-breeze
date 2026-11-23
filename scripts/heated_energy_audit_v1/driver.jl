include("../omega3_long_v1/model.jl")
include("../omega3_long_v1/initialization.jl")
@inline external_heating(i,j,k,F,p)=p.heat*(conduction(i,j,k,F.T,p)-F.ρᵈ[i,j,k]*(F.T[i,j,k]-p.tt)*cool_rate(k,p))
include("audit.jl")
function main()
    out=ARGS[1];mkpath(out);root="data/omega3-long-v1-run-2096858/parent"
    cfg=TOML.parsefile(joinpath(root,"metadata.toml"));cfg["initial_rms"]=0.;cfg["perturbation"]=0.
    small=merge(cfg,Dict("nx"=>8,"ny"=>8,"nz"=>24));m,p=build_model(small;arch=CPU(),perturb=false)
    set!(m;u=(x,y,z)->.03sinpi(x/2)*cospi(y/2),v=(x,y,z)->.02cospi(x/2)*sinpi(y/2),w=0.)
    Oceananigans.TimeSteppers.update_state!(m);a=audit(m,p)
    @assert abs(a["mass_tendency"])<1e-10 && abs(a["entropy_heat"])<1e-10 && abs(a["thermal_source_identity_error"])<1e-10
    write_toml(joinpath(out,"CHECKS.toml"),Dict("status"=>"PASS","periodic_mass_tendency"=>a["mass_tendency"],"closed_entropy_heat"=>a["entropy_heat"],"thermal_identity"=>a["thermal_source_identity_error"]))
    m,p=build_model(cfg;arch=CPU(),perturb=false)
    file=joinpath(root,"fields3d/state_0000800.jld2");x=load(file);rho=Float64.(x["density"]);T=Float64.(x["temperature"])
    set!(m;ρ=rho,θ=(rho.*T).^.6 ./rho,u=Float64.(x["u_face"]),v=Float64.(x["v_face"]),w=Float64.(x["w_face"]),qᵗ=0.)
    Oceananigans.TimeSteppers.update_state!(m);a=audit(m,p)
    lines=readlines(joinpath(root,"diagnostics.csv"));keys=split(lines[1],',');rr=[parse.(Float64,split(l,',')) for l in lines[2:end]]
    ti=findfirst(==("time"),keys);ei=findfirst(==("total_energy"),keys);idx=argmin([abs(r[ti]-800) for r in rr]);fd=(rr[idx+1][ei]-rr[idx-1][ei])/(rr[idx+1][ti]-rr[idx-1][ti])
    a["saved_time"]=x["t"];a["recorded_energy_derivative"]=fd;a["relative_derivative_mismatch"]=abs(fd-a["predicted_energy_derivative"])/max(abs(fd),1e-20)
    a["bottom_heat_input"]=16p.flux;a["defect_over_bottom_input"]=a["semidiscrete_energy_defect"]/(16p.flux)
    write_toml(joinpath(out,"HEATED_AUDIT.toml"),a);println(a)
end
main()

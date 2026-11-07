include("run.jl")
include("check_output.jl")
using Test

function check_diagnostics(root; heating=false)
    lines=readlines(joinpath(root,"diagnostics.csv"))
    names=split(first(lines), ',')
    data=permutedims(reduce(hcat, [parse.(Float64,split(line, ',')) for line in lines[2:end]]))
    col(name)=data[:,findfirst(==(name),names)]
    @test all(isfinite,data)
    @test all(diff(col("time_seconds")) .> 0)
    @test all(col("kinetic_energy_J_m2") .>= 0)
    @test all(col("sponge_power_W_m2") .<= 1e-12)
    if heating
        measured=last(col("theta_column_kg_K_m2"))-first(col("theta_column_kg_K_m2"))
        expected=last(col("bottom_theta_flux_kg_K_m2_s"))*(last(col("time_seconds"))-first(col("time_seconds")))
        @test isapprox(measured,expected;rtol=1e-4)
    end
    println("Diagnostic checks passed: ",root)
end

root=length(ARGS)>0 ? ARGS[1] : "data/convection48-v1"
arch=length(ARGS)>1 ? ARGS[2] : "GPU"
# New diagnostics must pass analytic rest/heating controls before the long runs.
for name in ("rest","heating")
    out=joinpath(root,"check_"*name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_diagnostics(out;heating=name=="heating")
    GC.gc()
end
for tag in ("p020","p002"), forcing in ("on","off")
    name="convection48_"*forcing*"_"*tag
    out=joinpath(root,name)
    main(["configs/"*name*".toml",out,arch])
    check_output(out)
    check_diagnostics(out)
    GC.gc()
end

using Pkg, TOML
const REV = "8273be89d9357ca2ac572347c45db8c647e3c45c"
project = ENV["TE_PROJECT"]
@assert !isfile(joinpath(project, "SETUP_PASS.toml"))
mkpath(project)
open(joinpath(project,"Project.toml"),"w") do io
    print(io, """
[deps]
Breeze = "660aa2fb-d4c8-4359-a52c-9c057bc511da"
Oceananigans = "9e8cae18-63c1-5223-a75c-80ca9d6e9a09"
CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
CUDA_Runtime_jll = "76a88914-d11a-5bdc-97e0-2f5a05c973a2"
CUDA_Compiler_jll = "d1e2174e-dfdc-576e-b43e-73b79eb1aca8"
[compat]
CUDA = "=6.4.2"
JLD2 = "=0.6.7"
Oceananigans = "=0.118.1"
julia = "1.12"
[preferences.CUDA_Runtime_jll]
version = "12.2"
[preferences.CUDA_Compiler_jll]
version = "12.2"
""")
end
Pkg.activate(project)
Pkg.add(PackageSpec(url="https://github.com/NumericalEarth/Breeze.jl", rev=REV))
Pkg.instantiate()
m=TOML.parsefile(joinpath(project,"Manifest.toml"))
@assert only(m["deps"]["Breeze"])["repo-rev"]==REV
open(joinpath(project,"SETUP_PASS.toml"),"w") do io
    TOML.print(io, Dict("status"=>"PASS", "revision"=>REV, "project"=>project))
end
println("ISOLATED_ENVIRONMENT_READY");flush(stdout)

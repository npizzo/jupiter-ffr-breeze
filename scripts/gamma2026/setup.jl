# Run on Unity. This isolated environment never resolves the Breeze environment.
using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "environments", "gamma2026"))
Pkg.add(PackageSpec(url="https://github.com/FourierFlows/GeophysicalFlows.jl", rev="51ca98a0d471ce4ef8bbba8b73daf2bad54d53fb"))
Pkg.instantiate()
Pkg.precompile()
Pkg.status()

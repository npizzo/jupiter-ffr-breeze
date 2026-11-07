using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.add(PackageSpec(url="https://github.com/NumericalEarth/Breeze.jl", rev="0479d9c8ce006adcb1e3314c324b84980ad4871b"))
Pkg.instantiate()
Pkg.precompile()

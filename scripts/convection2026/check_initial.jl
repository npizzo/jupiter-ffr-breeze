using Test, Statistics
include("common.jl")
include("initial_flow.jl")
root=only(ARGS);mkpath(root)
@testset "New resolution campaign initialization" begin
    for n in (512,1024)
        reference=nothing
        for heat in HEATS
            c=TOML.parsefile(config_path(heat));c["nx"]=n;c["ny"]=n;c["nz"]=2
            @test c["jet_speed"]==0 && c["shear_rossby"]==0 && c["beta"]==0
            @test c["latitude"]==85 && c["heat_flux"]==heat
            f=2c["rotation_rate"]*sind(c["latitude"])
            u,v,w,p,info=initial_large_eddies(c,[1.0,0.5],f)
            @test all(isfinite,u) && all(isfinite,v) && all(isfinite,p)
            @test maximum(abs,w)==0 && maximum(abs,v[:,[1,end],:])==0
            @test info["initial_mass_divergence_relative"]<1e-12
            @test info["initial_random_geostrophic_residual_f0"]<1e-12
            @test u[:,:,1]==u[:,:,2] && v[:,:,1]==v[:,:,2]
            speed=sqrt(mean(u.^2 .+(v[:,1:end-1,:].^2 .+v[:,2:end,:].^2)./2))
            @test isapprox(speed/(abs(f)*c["turbulence_wavelength"]/(2pi)),0.1;rtol=1e-12)
            if reference!==nothing; @test (u,v,w,p)==reference;end
            reference=(u,v,w,p)
        end
    end
    for f in readdir(@__DIR__;join=true)
        endswith(f,".jl") && Meta.parseall(read(f,String))
    end
end
run(`tools/field-export-env/bin/python scripts/convection2026/transfer.py --self-test`)
write_gate(joinpath(root,"UNIT_PASS.toml"))

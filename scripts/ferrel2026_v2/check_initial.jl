using Test, Statistics
include("common.jl")
include("initial_flow.jl")
include("jet_initialization.jl")
root=only(ARGS);mkpath(root)
@testset "Midlatitude configuration and initialization" begin
    for n in (96,512,1024)
        reference=nothing
        for case in ("b0","b162","boundary162")
            c=TOML.parsefile(config_path(case));c["nx"]=n;c["ny"]=n;c["nz"]=2
            n==96 && (c["turbulence_wavelength"]=3200000.0)
            @test c["jet_speed"]==0 && c["shear_rossby"]==0
            @test c["latitude"]==55 && c["rotation_model"]=="beta-plane"
            @test c["beta"]≈2c["rotation_rate"]*cosd(c["latitude"])/c["planet_radius"]
            f=2c["rotation_rate"]*sind(c["latitude"])
            u,v,w,p,info=initial_large_eddies(c,[1.0,0.5],f)
            jet,pressure,jet_info=balanced_cosine_jet(c,f)
            @test jet_info["initial_jet_curvature_over_beta"]<1
            @test jet_info["initial_jet_pressure_residual"]<1e-15
            @test abs(mean(jet))<1e-12
            @test maximum(abs,w)==0 && maximum(abs,v[:,[1,end],:])==0
            @test info["initial_mass_divergence_relative"]<1e-12
            @test info["initial_random_geostrophic_residual_f0"]<1e-12
            speed=sqrt(mean(u.^2 .+(v[:,1:end-1,:].^2 .+v[:,2:end,:].^2)./2))
            @test isapprox(speed/(abs(f)*c["turbulence_wavelength"]/(2pi)),0.1;rtol=1e-12)
            if reference!==nothing;@test (u,v,w,p,jet,pressure)==reference;end
            reference=(u,v,w,p,jet,pressure)
        end
    end
    for case in ("axis0","axis162")
        c=TOML.parsefile(config_path(case))
        @test c["perturbation"]==0 && c["initial_velocity"]=="rest"
    end
    for f in readdir(@__DIR__;join=true)
        endswith(f,".jl") && Meta.parseall(read(f,String))
    end
end
run(`tools/field-export-env/bin/python scripts/ferrel2026_v2/transfer.py --self-test`)
write_gate(joinpath(root,"UNIT_PASS.toml"))

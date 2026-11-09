using Test, TOML, SHA
include("initial_flow.jl")
root=ARGS[1]; mkpath(root)
cfg=TOML.parsefile("configs/eddy20_f_heated.toml")
merge!(cfg,Dict("nx"=>128,"ny"=>128,"nz"=>8))
f=2cfg["rotation_rate"]*sind(cfg["latitude"])
rho=exp.(-collect(1:8)./8)
@testset "Large random eddy constraints" begin
    u,v,w,p,info=initial_large_eddies(cfg,rho,f)
    @test all(isfinite,u) && all(isfinite,v) && all(isfinite,p)
    @test maximum(abs,w)==0
    @test maximum(abs,v[:,[1,end],:])==0
    @test info["initial_mass_divergence_relative"]<1e-12
    @test info["initial_random_geostrophic_residual_f0"]<1e-12
    @test u[:,:,1]==u[:,:,end] && v[:,:,1]==v[:,:,end]
    @test maximum(abs,mean(u,dims=1))<1e-12
    @test maximum(abs,mean(v,dims=1))<1e-12
    rms=sqrt(mean(u.^2 .+(v[:,1:end-1,:].^2 .+v[:,2:end,:].^2)./2))
    @test isapprox(rms/(abs(f)*cfg["turbulence_wavelength"]/(2pi)),0.1;rtol=1e-13)
    @test info["initial_seed_wavelength_min_m"] >= cfg["turbulence_wavelength"]/1.25
    @test info["initial_seed_wavelength_max_m"] <= cfg["turbulence_wavelength"]/0.75
    for plane in ("f","beta"), heating in ("heated","control")
        c=TOML.parsefile("configs/eddy20_$(plane)_$(heating).toml")
        merge!(c,Dict("nx"=>128,"ny"=>128,"nz"=>8))
        uu,vv,ww,pp,_=initial_large_eddies(c,rho,f)
        @test (u,v,w,p)==(uu,vv,ww,pp)
    end
    c=copy(cfg);c["turbulence_seed"]+=1
    uu,_,_,_,_=initial_large_eddies(c,rho,f)
    @test u!=uu
    # Same physical modes at production resolution, still cheap with two levels.
    c=copy(cfg);merge!(c,Dict("nx"=>1024,"ny"=>1024,"nz"=>2))
    _,_,_,_,fine=initial_large_eddies(c,[1.0,0.8],f)
    @test fine["initial_horizontal_mode_count"]==info["initial_horizontal_mode_count"]
    @test fine["initial_mass_divergence_relative"]<1e-12
    record=merge(info,Dict("status"=>"PASS","initializer_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"initial_flow.jl"))))))
    open(joinpath(root,"UNIT_PASS.toml"),"w") do io; TOML.print(io,record); end
end
for file in readdir(@__DIR__;join=true)
    endswith(file,".jl") && Meta.parseall(read(file,String))
end

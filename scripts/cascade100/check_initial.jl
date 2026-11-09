using Test,TOML,SHA
include("initial_flow.jl")
root=only(ARGS);mkpath(root)
cfg=TOML.parsefile("configs/cascade100_f_heated.toml")
cfg["nz"]=2;rho=[1.0,0.5];f=2cfg["rotation_rate"]*sind(cfg["latitude"])
@testset "Coarse no-shear initial conditions" begin
    @test cfg["nx"]==256 && cfg["ny"]==256
    @test cfg["shear_rossby"]==0 && cfg["beta"]==0 && cfg["jet_speed"]==0
    u,v,w,p,info=initial_large_eddies(cfg,rho,f)
    @test all(isfinite,u) && all(isfinite,v) && all(isfinite,p)
    @test maximum(abs,w)==0
    @test maximum(abs,v[:,[1,end],:])==0
    @test info["initial_mass_divergence_relative"]<1e-12
    @test info["initial_random_geostrophic_residual_f0"]<1e-12
    @test u[:,:,1]==u[:,:,end] && v[:,:,1]==v[:,:,end]
    rms=sqrt(mean(u.^2 .+(v[:,1:end-1,:].^2 .+v[:,2:end,:].^2)./2))
    @test isapprox(rms/(abs(f)*cfg["turbulence_wavelength"]/(2pi)),.1;rtol=1e-12)
    c=TOML.parsefile("configs/cascade100_f_control.toml");c["nz"]=2
    uu,vv,ww,pp,_=initial_large_eddies(c,rho,f)
    @test (u,v,w,p)==(uu,vv,ww,pp)
    for file in readdir(@__DIR__;join=true)
        endswith(file,".jl") && Meta.parseall(read(file,String))
    end
end
run(`tools/field-export-env/bin/python scripts/cascade100/transfer.py --self-test`)
open(joinpath(root,"UNIT_PASS.toml"),"w") do io
 TOML.print(io,Dict("status"=>"PASS","initializer_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"initial_flow.jl")))),"transfer_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"transfer.py"))))))
end

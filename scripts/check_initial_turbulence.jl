using Test, TOML

function check_initial_turbulence()
    cfg=TOML.parsefile("configs/turbulence_check.toml")
    f=2cfg["rotation_rate"]*sind(cfg["latitude"])
    rho=exp.(-collect(1:cfg["nz"])./cfg["nz"])
    u,v,w,info=initial_turbulence(cfg,rho,f)
    @test all(isfinite,u) && all(isfinite,v) && all(isfinite,w)
    @test info["initial_mass_divergence_relative"] < 1e-12
    @test isapprox(info["initial_eddy_rossby"],0.1;rtol=1e-12)
    @test isapprox(info["initial_vertical_rms_lower_m_s"]/info["initial_horizontal_rms_lower_m_s"],0.2;rtol=1e-12)
    @test maximum(abs,w[:,:,[1,end]]) == 0
    u2,v2,w2,_=initial_turbulence(cfg,rho,f)
    @test u==u2 && v==v2 && w==w2
    cfg["turbulence_seed"]+=1
    u3,_,_,_=initial_turbulence(cfg,rho,f)
    @test u != u3
    println("Initial turbulence constraint, amplitude and reproducibility checks passed")
end

function check_turbulent_output(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    initial=load(joinpath(out,"frame_00000.jld2"))
    lower=initial["z"] .< m["cap_height"]
    h_rms=sqrt(mean(2 .* initial["horizontal_ke_profile"][lower]./initial["density"][lower]))
    w_rms=sqrt(mean(2 .* initial["vertical_ke_profile"][lower]./initial["density"][lower]))
    @test isapprox(h_rms,m["initial_horizontal_rms_lower_m_s"];rtol=1e-9)
    @test isapprox(w_rms,m["initial_vertical_rms_lower_m_s"];rtol=1e-9)
    @test m["initial_mass_divergence_relative"] < 1e-12
    println("Initial turbulence in solver output matches prescribed amplitudes: ",out)
end

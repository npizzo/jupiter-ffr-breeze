using Test, TOML, JLD2, Statistics
include("../ferrel2026_v3/check_ferrel.jl")

function check_seed_case(out,cfg,reference_path)
    check_ferrel_output(out)
    metadata=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test metadata["status"]=="COMPLETE"
    @test metadata["simulated_seconds"]==21600.0
    @test metadata["segment_start_seconds"]==0.0
    @test metadata["pickup_checkpoint"]==""
    @test metadata["thermal_seed_record_sha256"]==seed_sha(joinpath(out,"thermal_seed.toml"))
    info=TOML.parsefile(joinpath(out,"thermal_seed.toml"))
    reference=TOML.parsefile(reference_path)
    @test info["model"]==cfg["thermal_seed_model"]
    @test info["reference_grid"]==[512,512,96]
    for key in ("reference_grid","reference_height_m","reference_seed","reference_raw_variance",
                "reference_target_variance_K2","amplitude_parameter_K","envelope",
                "phase_denominator","modes")
        @test info[key]==reference[key]
    end
    if cfg["thermal_seed_model"]=="zero"
        @test all(info["actual_preaddition_variance_K2"].==0)
    else
        @test all(isapprox.(info["actual_preaddition_variance_K2"],info["reference_target_variance_K2"];rtol=1024eps(Float64),atol=0))
    end
    planes=metadata["precision_diagnostic_indices_julia_1based"]
    frames=sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))
    @test length(frames)==7
    for (index,path) in enumerate(frames)
        jldopen(path) do d
            @test d["t"]==(index-1)*3600.0
            @test d["precision_indices_julia_1based"]==planes
            @test d["precision_heights"]==metadata["precision_diagnostic_heights_actual_m"]
            w=d["precision_wcenter_slices"];theta=d["precision_theta_anomaly_slices"]
            @test eltype(w)==eltype(theta)==Float64
            @test size(w)==size(theta)==(cfg["nx"],cfg["ny"],5)
            @test all(isfinite,w) && all(isfinite,theta)
            theta_profile=d["theta_profile"];rho=d["density"];pressure=d["pressure"]
            for (q,k) in enumerate(planes)
                # Derived reduction bounds, not a relative tolerance against a
                # near-zero heat flux. The numerator bounds Float64 arithmetic
                # in two independent means and the same-plane subtraction.
                n=cfg["nx"]*cfg["ny"]
                gamma=(4n+32)*eps(Float64)/(1-(4n+32)*eps(Float64))
                power=mean(abs2,w[:,:,q])
                expected_power=2d["vertical_ke_profile"][k]/rho[k]
                @test abs(power-expected_power) <= gamma*(abs(power)+abs(expected_power)) + floatmin(Float64)
                covariance=mean(w[:,:,q].*theta[:,:,q])
                exner=(pressure[k]/1e5)^((8.314462618/cfg["molar_mass"])/cfg["heat_capacity"])
                expected_covariance=d["heat_flux_profile"][k]/(rho[k]*cfg["heat_capacity"]*exner)
                allowance=gamma*(mean(abs,w[:,:,q].*theta[:,:,q])+abs(expected_covariance))+floatmin(Float64)
                @test abs(covariance-expected_covariance) <= allowance
                v=findfirst(==(d["precision_heights"][q]),d["video_heights"])
                if v!==nothing
                    @test Float32.(w[:,:,q])==d["wcenter_slices"][:,:,v]
                end
                if index==1
                    expected=cfg["thermal_seed_model"]=="zero" ? zeros(Float64,cfg["nx"],cfg["ny"]) : horizontal_seed_plane(cfg,k,info["reference_target_variance_K2"][k])
                    # The model initializes absolute Float64 theta through
                    # density/thermodynamic conversions.16 ulps is an explicit
                    # initialization conversion allowance, not weak-field
                    # transport accuracy or a later-time spectral tolerance.
                    allowance=16eps(theta_profile[k])
                    @test maximum(abs,theta[:,:,q].-expected) <= allowance
                end
            end
        end
    end
    println("PASS: fresh6h seed case, native v3 physics checks, exact precision-plane indices and Float64 transport consistency")
end

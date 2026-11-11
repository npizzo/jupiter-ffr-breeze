include("check_output.jl")
include("restart_validation.jl")
function check_ferrel_output(out)
    check_shear_output(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"))
    @test m["native_budget_checks"]>0
    @test m["initial_jet_curvature_over_beta"]<1
    paths=sort(filter(p->startswith(basename(p),"frame_")&&endswith(p,".jld2"),readdir(out;join=true)))
    for path in paths
        jldopen(path) do d
            @test size(d["native_rhov"])==(m["ny"]+1,m["nz"])
            @test size(d["native_rhow"])==(m["ny"],m["nz"]+1)
            @test size(d["mass_streamfunction"])==(m["ny"]+1,m["nz"]+1)
            @test size(d["native_rhou_rhs"])==(m["ny"],m["nz"],4)
            for key in ("native_rhov","native_rhow","mass_streamfunction","native_rhou_rhs","centered_native_advection_remainder","eddy_to_mean_ke_transfer_y","eddy_to_mean_ke_transfer_z")
                @test all(isfinite,d[key])
            end
            massscale=max(maximum(abs,d["native_rhov"]),maximum(abs,d["native_rhow"]),1e-10)
            @test d["mass_flux_reconstruction_max"] <= 1e-6massscale+1e-9
            @test d["mass_divergence_max"] <= 1e-6massscale/(m["height"]/m["nz"])+1e-12
            @test maximum(abs,d["mass_streamfunction"][:,1])==0
            @test maximum(abs,d["native_rhov"][[1,end],:])<1e-10
            @test maximum(abs,d["native_rhow"][:,[1,end]])<1e-10
            @test d["total_rho_vtracer"]≈d["mean_rho_vtracer"]+d["zonal_rho_vtracer"]
            @test d["total_rho_wtracer"]≈d["mean_rho_wtracer"]+d["zonal_rho_wtracer"]
            if startswith(m["name"],"ferrel2026_axis")
                # Exact x symmetry is a control, including the thermal initial state.
                for key in ("u","v","w","theta")
                    a=d[key]
                    @test maximum(abs,a.-mean(a;dims=1))<1e-8
                end
                @test abs(d["eddy_to_mean_ke_transfer_W_m2"])<1e-8
            end
        end
    end
    rows=readlines(joinpath(out,"native_budget","checks.csv"))
    @test length(rows)==m["native_budget_checks"]+1
    println("Native momentum closure, mass circulation, and tracer decomposition passed: ",out)
end

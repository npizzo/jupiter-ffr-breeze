using Test, Statistics, TOML
include("common.jl")
include("thermal_seed.jl")
include("../ferrel2026_v3/common.jl")

require_allocation()
length(ARGS)==2 || error("usage: validate_initialization.jl FRESH_GATE_ROOT V3_PREFLIGHT_ROOT")
root,prior=ARGS
verify_gate(joinpath(prior,"PASS512.toml"))
fresh_directory(root)
audit=audit_seed_source()
write(joinpath(root,"source_inheritance_audit.json"),audit)
cfg=TOML.parsefile(seed_config_path("modes512_b162"))
target,raw=reference_seed_variances(cfg)
H,Lx,Ly=cfg["height"],cfg["lx"],cfg["ly"]
Nz=cfg["nz"]; z=((1:Nz).-0.5).*H/Nz
analytic_roundoff=1024eps(Float64)
worst_variance_relative=0.0
worst_normalized_mean=0.0
worst_normalized_point_error=0.0

@testset "Analytic horizontal seed and original variance reference" begin
    @test all(target .> 0)
    @test target == cfg["perturbation"]^2 .* seed_envelope.(z,H).^2 .* raw
    @test cfg["perturbation"]==1e-5
    @test all(375000.0 <= 2pi/hypot(2pi*q.m/Lx,pi*q.n/Ly) <= 750000.0 for q in THERMAL_MODES)
    for k in (1,8,48,96),x in (0.0,0.137Lx,0.73Lx),y in (-Ly/2,-0.11Ly,Ly/2)
        @test abs(modal_value(x+Lx,y,k,Lx,Ly)-modal_value(x,y,k,Lx,Ly)) <= analytic_roundoff*sum(abs(q.a) for q in THERMAL_MODES)
        for wall in (-Ly/2,Ly/2)
            derivative=sum(-q.a*pi*q.n/Ly*cos(2pi*q.m*x/Lx+phase_radians(q,k))*sin(pi*q.n*(wall+Ly/2)/Ly) for q in THERMAL_MODES)
            @test abs(derivative) <= analytic_roundoff*sum(abs(q.a*pi*q.n/Ly) for q in THERMAL_MODES)
        end
    end
    #64/128 are analytic sampling/orthogonality checks only: no model evolution,
    # domain change, wavelength adaptation or physical-resolution assertion.
    #512/1024 test every production z-plane at the actual physical wavelengths.
    for n in (64,128,512,1024),k=1:Nz
        a=horizontal_seed_plane(cfg,k,target[k];nx=n,ny=n)
        variance_error=abs(mean(abs2,a)/target[k]-1)
        mean_error=abs(mean(a))/sqrt(target[k])
        @test variance_error <= analytic_roundoff
        @test mean_error <= analytic_roundoff
        global worst_variance_relative=max(worst_variance_relative,variance_error)
        global worst_normalized_mean=max(worst_normalized_mean,mean_error)
        # Coarse and fine centers differ. Test each against the same continuous
        # formula rather than falsely equating every second fine cell to coarse.
        for (i,j) in ((1,1),(n÷3,2n÷3),(n÷2,n÷2),(n,n))
            x=(i-0.5)*Lx/n;y=(j-0.5)*Ly/n-Ly/2
            expected=sqrt(target[k]/modal_variance())*modal_value(x,y,k,Lx,Ly)
            discrepancy=abs(a[i,j]-expected)/sqrt(target[k])
            @test discrepancy <= analytic_roundoff
            global worst_normalized_point_error=max(worst_normalized_point_error,discrepancy)
        end
    end
    small=merge(cfg,Dict("nx"=>64,"ny"=>64))
    theta_ref(z)=cfg["theta0"]*exp((cfg["convective_n2"]*min(z,cfg["cap_height"])+cfg["cap_n2"]*max(z-cfg["cap_height"],0))/cfg["gravity"])
    theta,info=build_seeded_theta(small,theta_ref)
    @test info["reference_target_variance_K2"]==target
    @test all(isapprox.(info["actual_preaddition_variance_K2"],target;rtol=analytic_roundoff,atol=0))
    open(joinpath(root,"thermal_mode_reference.toml"),"w") do io;TOML.print(io,info);end
    cached=merge(small,Dict("thermal_seed_reference_file"=>abspath(joinpath(root,"thermal_mode_reference.toml")),
        "thermal_seed_reference_sha256"=>seed_sha(joinpath(root,"thermal_mode_reference.toml"))))
    @test reference_seed_variances(cached)==(target,raw)
    @test_throws ErrorException reference_seed_variances(merge(cached,Dict("thermal_seed_reference_sha256"=>"invalid")))
    zero_cfg=merge(small,Dict("thermal_seed_model"=>"zero"))
    zero_theta,zero_info=build_seeded_theta(zero_cfg,theta_ref)
    @test all(zero_info["actual_preaddition_variance_K2"].==0)
    for k=1:Nz
        @test all(zero_theta[:,:,k].==theta_ref(z[k]))
        # Float64 absolute-temperature addition only. This tests construction,
        # not subsequent liquid-theta conversions in the model.
        original_plane=horizontal_seed_plane(small,k,target[k])
        @test maximum(abs,(theta[:,:,k].-theta_ref(z[k])).-original_plane) <= eps(theta_ref(z[k]))
    end
    for path in readdir(@__DIR__;join=true)
        endswith(path,".jl") && Meta.parseall(read(path,String))
    end
end
write_seed_gate(joinpath(root,"INITIALIZATION_PASS.toml");
    prior_pass512_sha256=seed_sha(joinpath(prior,"PASS512.toml")),
    runtime_marker_sha256=seed_sha(joinpath(prior,"runtime.toml")),
    analytic_grids=[64,128],production_initialization_grids=[512,1024],z_planes_tested=Nz,
    max_relative_variance_error=worst_variance_relative,max_mean_over_rms=worst_normalized_mean,
    max_point_error_over_rms=worst_normalized_point_error,
    thermal_mode_reference_sha256=seed_sha(joinpath(root,"thermal_mode_reference.toml")),
    interpretation="Initialization/source-inheritance validation only. No new dynamical run, convergence or convective-fidelity claim; prior v3 physics checks remain separately identified.")
println("PASS: fixed physical horizontal modes, both full-grid initializations, exact original512 variance reference, zero control, and source inheritance")

using Random, Statistics, TOML, SHA

# x-periodic / y-Neumann analytic modes. All constants are part of provenance.
# The coefficient and phase table is independent of Nx/Ny. Layer-dependent
# phases deliberately retain short vertical correlations: this is a horizontal
# seed test, not a vertically smooth or vertically converged initialization.
const THERMAL_MODES = [
    (m=16, n=0,  a=1.0,  phase0=17,  stride=431),
    (m=20, n=8,  a=-1.0, phase0=131, stride=613),
    (m=24, n=16, a=1.0,  phase0=281, stride=337),
    (m=28, n=16, a=-1.0, phase0=419, stride=719),
    (m=16, n=32, a=1.0,  phase0=557, stride=227),
    (m=24, n=32, a=-1.0, phase0=683, stride=829),
    (m=8,  n=40, a=1.0,  phase0=797, stride=491),
    (m=12, n=48, a=-1.0, phase0=919, stride=157),
]
const PHASE_DENOMINATOR = 997
phase_numerator(mode,k) = mod(mode.phase0 + mode.stride*(k-1),PHASE_DENOMINATOR)
phase_radians(mode,k) = 2pi*phase_numerator(mode,k)/PHASE_DENOMINATOR
modal_variance() = sum(q.a^2*(q.n==0 ? 0.5 : 0.25) for q in THERMAL_MODES)
seed_envelope(z,H) = exp(-(z/(0.35H))^2)

"""Replay the original512 cellwise seed exactly before adding theta_ref.

The reference draw and horizontal demeaning use the same Float64 array shape,
RNG, order and operations as frozen v3. Return realized variance, not A²/12.
This reference is always512²×Nz, including when the target grid is1024².
"""
function reference_seed_variances(cfg)
    cfg["thermal_seed_reference_nx"] == 512 || error("Reference grid must remain512")
    cfg["thermal_seed_reference_ny"] == 512 || error("Reference grid must remain512")
    cfg["seed"] == 42 || error("Reference seed must remain42")
    H,Nz=cfg["height"],cfg["nz"]
    if haskey(cfg,"thermal_seed_reference_file")
        path=cfg["thermal_seed_reference_file"]
        bytes2hex(sha256(read(path)))==cfg["thermal_seed_reference_sha256"] || error("Changed recorded variance reference")
        reference=TOML.parsefile(path)
        reference["reference_grid"]==[512,512,Nz] || error("Variance reference grid mismatch")
        reference["reference_height_m"]==H || error("Variance reference height mismatch")
        reference["reference_seed"]==cfg["seed"] || error("Variance reference RNG mismatch")
        reference["amplitude_parameter_K"]==cfg["perturbation"] || error("Variance reference amplitude mismatch")
        # Production consumes the exact CPU-gated table. It does not recompute
        # reductions on another node and silently perturb the matched RMS.
        return Float64.(reference["reference_target_variance_K2"]),Float64.(reference["reference_raw_variance"])
    end
    noise=rand(MersenneTwister(cfg["seed"]),Float64,512,512,Nz).-0.5
    noise .-= mean(noise,dims=(1,2))
    raw_variance=vec(mean(abs2,noise;dims=(1,2)))
    noise=nothing
    z=((1:Nz).-0.5).*H/Nz
    target=cfg["perturbation"]^2 .* seed_envelope.(z,H).^2 .* raw_variance
    return target,raw_variance
end

function modal_value(x,y,k,Lx,Ly)
    sum(q.a*cos(2pi*q.m*x/Lx+phase_radians(q,k))*
        cos(pi*q.n*(y+Ly/2)/Ly) for q in THERMAL_MODES)
end

function horizontal_seed_plane(cfg,k,target_variance;nx=cfg["nx"],ny=cfg["ny"])
    Lx,Ly=cfg["lx"],cfg["ly"]
    maximum(q.m for q in THERMAL_MODES) < nx/2 || error("Zonal seed reaches Nyquist")
    maximum(q.n for q in THERMAL_MODES) < ny || error("Meridional seed reaches Nyquist")
    x=((1:nx).-0.5).*Lx/nx
    y=((1:ny).-0.5).*Ly/ny.-Ly/2
    scale=sqrt(target_variance/modal_variance())
    a=zeros(Float64,nx,ny)
    for q in THERMAL_MODES
        cx=cos.(2pi*q.m.*x./Lx .+ phase_radians(q,k))
        cy=cos.(pi*q.n.*(y.+Ly/2)./Ly)
        a .+= (scale*q.a) .* reshape(cx,nx,1) .* reshape(cy,1,ny)
    end
    return a
end

function build_seeded_theta(cfg,theta_ref)
    model=cfg["thermal_seed_model"]
    model in ("zero","horizontal_modes") || error("Unknown thermal seed")
    Nx,Ny,Nz=cfg["nx"],cfg["ny"],cfg["nz"]
    H=cfg["height"]; z=((1:Nz).-0.5).*H/Nz
    target,raw=reference_seed_variances(cfg)
    theta=Array{Float64}(undef,Nx,Ny,Nz)
    actual=Float64[]; plane_mean=Float64[]
    for k=1:Nz
        plane=model=="zero" ? zeros(Float64,Nx,Ny) : horizontal_seed_plane(cfg,k,target[k])
        push!(actual,mean(abs2,plane));push!(plane_mean,mean(plane))
        theta[:,:,k] .= theta_ref(z[k]) .+ plane
    end
    info=Dict{String,Any}(
        "model"=>model,"reference_grid"=>[512,512,Nz],"reference_height_m"=>H,"reference_seed"=>cfg["seed"],
        "reference_raw_variance"=>raw,"reference_target_variance_K2"=>target,
        "actual_preaddition_variance_K2"=>actual,"actual_preaddition_mean_K"=>plane_mean,
        "amplitude_parameter_K"=>cfg["perturbation"],"envelope"=>"exp(-(z/(0.35H))^2)",
        "variance_scope"=>"Exact realized original512 seed42 variance after horizontal demeaning, amplitude and envelope, before adding absolute theta_ref; zero control deliberately omits it",
        "horizontal_grid_independent"=>true,"vertical_smoothness"=>"none claimed; deterministic layer-dependent phases at unchanged Nz96",
        "phase_denominator"=>PHASE_DENOMINATOR,
        "modes"=>[Dict("m"=>q.m,"n"=>q.n,"coefficient"=>q.a,"phase0_numerator"=>q.phase0,
                        "layer_stride_numerator"=>q.stride,
                        "phase_numerators_by_layer"=>[phase_numerator(q,k) for k=1:Nz],
                        "wavelength_m"=>2pi/hypot(2pi*q.m/cfg["lx"],pi*q.n/cfg["ly"])) for q in THERMAL_MODES])
    return theta,info
end

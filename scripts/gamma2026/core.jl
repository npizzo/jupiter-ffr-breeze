module Gamma2026
using GeophysicalFlows, CUDA, FFTW, Random, Statistics, LinearAlgebra, TOML, JLD2, SHA, Printf
const QG = GeophysicalFlows.SingleLayerQG
const GF_COMMIT = "51ca98a0d471ce4ef8bbba8b73daf2bad54d53fb"
const DAY = 86400.0

function scales(c)
    U=Float64(c["rms_velocity_m_s"]); gamma=Float64(c["gamma_m2_s"])
    L=(U/gamma)^(1/3)
    return (;U,gamma,L,T=L/U)
end
source_hash()=bytes2hex(sha256(read(@__FILE__)))
config_hash(c)=bytes2hex(sha256(sprint(io->TOML.print(io,c;sorted=true))))

function problem(c, device; initial=true)
    n=Int(c["n"]); s=scales(c); ld=Float64(c["domain_Lgamma"])
    radius=Float64(c["trap_Lgamma"]); kind=c["kind"]
    kind in ("polar","flat","none") || error("Unknown PV geometry")
    eta(x,y)= x*x+y*y < radius^2 ? (kind=="polar" ? -(x*x+y*y)/2 :
                     kind=="flat" ? -radius^2/2 : 0.0) : 0.0
    p=QG.Problem(device; nx=n,Lx=ld,eta,β=0.0,deformation_radius=Inf,U=0.0,
        ν=0.0,μ=0.0,dt=1e-4,stepper=get(c,"stepper","FilteredRK4"),
        aliased_fraction=Float64(get(c,"aliased_fraction",1/3)),T=Float64)
    initial && initialize!(p,c)
    return p
end

function initialize!(p,c)
    n=p.grid.nx
    kh=sqrt.(Array(p.grid.Krsq)); s=scales(c)
    k0=s.L/Float64(c["initial_inverse_wavenumber_m"])
    halfwidth=Float64(c["annulus_fraction"])*k0
    k0+halfwidth < (n/3)*2pi/p.grid.Lx || error("Initial band is not dealiased")
    rng=MersenneTwister(Int(c["seed"]))
    # Transform real noise so Hermitian symmetry is exact, including kx=0.
    qh=FFTW.rfft(randn(rng,n,n))
    qh .*= abs.(kh.-k0).<=halfwidth
    q=FFTW.irfft(qh,n)
    x=Array(p.grid.x); y=Array(p.grid.y)
    r=sqrt.(reshape(x,:,1).^2 .+ reshape(y,1,:).^2)
    a=Float64(c["taper_inner_Lgamma"]); b=Float64(c["taper_outer_Lgamma"])
    0 <= a < b < c["trap_Lgamma"] || error("Taper must lie inside trap")
    taper=ifelse.(r.<=a,1.0,ifelse.(r.>=b,0.0,0.5.*(1 .+ cos.(pi.*(r.-a)./(b-a)))))
    q .*= taper
    q .-= mean(q)
    A=device_array(p.grid.device)
    QG.set_q!(p,A(q))
    e=QG.energy(p)
    e>0 || error("Empty initial spectrum")
    p.sol .*= sqrt(0.5/e) # Domain RMS speed is one nondimensionally.
    QG.updatevars!(p)
    return nothing
end

function metrics(p)
    QG.updatevars!(p)
    E=Float64(QG.energy(p)); Z=Float64(QG.enstrophy(p))
    QG.updatevars!(p)
    umax=Float64(maximum(abs,p.vars.u)); vmax=Float64(maximum(abs,p.vars.v))
    qmax=Float64(maximum(abs,p.vars.q))
    return (;E,Z,umax,vmax,qmax)
end

function safe_step!(p,c; target=Inf)
    m=metrics(p)
    dt=min(Float64(c["max_dt_seconds"])/scales(c).T,
        Float64(c["cfl"])*p.grid.dx/max(m.umax+m.vmax,1e-12),target-p.clock.t)
    isfinite(dt) && dt>0 || error("Invalid timestep")
    p.clock.dt=dt
    stepforward!(p)
    return nothing
end

function checkpoint!(p,c,path,initial_energy)
    tmp=path*".tmp"
    jldsave(tmp;sol=Array(p.sol),t=Float64(p.clock.t),step=p.clock.step,
        dt=Float64(p.clock.dt),initial_energy,config_hash=config_hash(c),source_hash=source_hash(),
        geophysicalflows_commit=GF_COMMIT)
    mv(tmp,path;force=true)
end

function restore!(p,c,path)
    d=JLD2.load(path)
    d["source_hash"]==source_hash() || error("Source changed since checkpoint")
    d["config_hash"]==config_hash(c) || error("Config changed since checkpoint")
    size(d["sol"])==size(p.sol) || error("Restart grid mismatch")
    p.sol .= device_array(p.grid.device)(d["sol"])
    p.clock.t=d["t"]; p.clock.step=d["step"];p.clock.dt=d["dt"]
    QG.updatevars!(p)
    return d["initial_energy"]
end

function snapshot!(p,c,path)
    QG.updatevars!(p)
    stride=max(1,div(p.grid.nx,Int(c["movie_grid"])))
    # Output reduction is publication-oriented, not a substituted dynamics grid.
    zeta=Float32.(Array(p.vars.q)[1:stride:end,1:stride:end])
    pv=zeta .+ Float32.(Array(p.params.eta)[1:stride:end,1:stride:end])
    s=scales(c)
    jldsave(path;time_seconds=p.clock.t*s.T,vorticity_s_inv=zeta.*Float32(s.U/s.L),
        pv_s_inv=pv.*Float32(s.U/s.L),x_m=Array(p.grid.x)[1:stride:end].*s.L,
        y_m=Array(p.grid.y)[1:stride:end].*s.L,dimension="2D barotropic; no physical depth samples",
        sampling_stride=stride)
end

function run(c,out; device=GPU(),resume=false,gate=nothing)
    if device isa GPU
        CUDA.functional() || error("CUDA unavailable: refusing CPU production fallback")
    end
    if gate!==nothing
        g=TOML.parsefile(gate)
        g["status"]=="PASSED" || error("GPU validation did not pass")
        g["source_hash"]==source_hash() || error("Gate belongs to changed source")
        g["n"]==c["n"] || error("Gate grid differs")
    end
    mkpath(out);snapdir=joinpath(out,"frames");mkpath(snapdir)
    metadata=merge(c,Dict("status"=>"RUNNING","model"=>"unforced barotropic QG",
       "reproduction_status"=>"RECONSTRUCTED_PILOT_NOT_EXACT_PAPER_REPLICATION",
       "limitations"=>"Original SI taper/grid/filter not verified; explicit reconstructed settings. 75-day pilot is not a mature-crystal test.",
       "source_hash"=>source_hash(),"config_hash"=>config_hash(c),
       "geophysicalflows_commit"=>GF_COMMIT,"compute_location"=>"Unity",
       "output_dimension"=>2,"depth_panels"=>"not applicable; no invented 3D depth copies"))
    metafile=joinpath(out,"metadata.toml")
    isfile(metafile) && !resume && error("Output already exists")
    p=problem(c,device);ck=joinpath(out,"checkpoint.jld2")
    E0=resume ? restore!(p,c,ck) : QG.energy(p)
    open(metafile,"w") do io;TOML.print(io,metadata);end
    s=scales(c); tend=Float64(c["end_days"])*DAY/s.T
    interval=Float64(c["output_days"])*DAY/s.T
    nextout=(floor(p.clock.t/interval+1e-8)+1)*interval
    diagfile=joinpath(out,"diagnostics.csv")
    if !resume
        open(diagfile,"w") do io
            println(io,"time_days,step,energy_m2_s2,pv_enstrophy_s2,umax_m_s,vmax_m_s,zeta_max_s_inv,energy_fraction,wall_seconds")
        end
        snapshot!(p,c,joinpath(snapdir,"frame_000000.jld2"))
    end
    wall0=time(); lastcheckpoint=wall0
    while p.clock.t < tend-1e-10
        safe_step!(p,c;target=min(tend,nextout))
        if p.clock.step%20==0 || p.clock.t>=nextout-1e-10
            m=metrics(p)
            all(isfinite,(m.E,m.Z,m.qmax)) || error("Nonfinite flow")
            m.E<=E0*Float64(c["max_energy_fraction"]) || error("Unforced energy growth exceeded tolerance")
            m.E>=E0*Float64(c["min_energy_fraction"]) || error("Excessive energy loss; resolution review required")
        end
        if p.clock.t>=nextout-1e-10 || p.clock.t>=tend-1e-10
            m=metrics(p); day=p.clock.t*s.T/DAY
            open(diagfile,"a") do io
                println(io,join((day,p.clock.step,m.E*s.U^2,m.Z*(s.U/s.L)^2,
                    m.umax*s.U,m.vmax*s.U,m.qmax*s.U/s.L,m.E/E0,time()-wall0),","))
            end
            frame=round(Int,p.clock.t/interval)
            snapshot!(p,c,joinpath(snapdir,@sprintf("frame_%06d.jld2",frame)))
            @printf("day=%.3f step=%d E/E0=%.7f wall=%.1fs\n",day,p.clock.step,m.E/E0,time()-wall0)
            flush(stdout);nextout+=interval
        end
        if time()-lastcheckpoint>300
            checkpoint!(p,c,ck,E0);lastcheckpoint=time()
        end
        if time()-wall0>Float64(c["max_wall_hours"])*3600
            checkpoint!(p,c,ck,E0)
            metadata["status"]="WALL_BUDGET_STOP";metadata["completed_days"]=p.clock.t*s.T/DAY
            open(metafile,"w") do io;TOML.print(io,metadata);end
            error("Wall budget reached with resumable checkpoint")
        end
    end
    checkpoint!(p,c,ck,E0)
    metadata["status"]="COMPLETE";metadata["completed_days"]=p.clock.t*s.T/DAY
    metadata["wall_seconds"]=time()-wall0
    metadata["energy_fraction"]=QG.energy(p)/E0
    open(metafile,"w") do io;TOML.print(io,metadata);end
    return metadata
end
end

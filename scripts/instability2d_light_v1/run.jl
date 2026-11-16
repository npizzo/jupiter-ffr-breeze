include("audit.jl")
using JLD2, SHA
out=ARGS[1];ispath(out)&&error("Refusing to overwrite existing output");mkpath(out)
seed=parse(Int,get(ENV,"EXPERIMENT_SEED","1701"))
target=parse(Float64,get(ENV,"TARGET_GROWTH_TIMES","120"))
n=128;dt=0.01;sigma=0.6;steps=round(Int,target/sigma/dt)
o=BandOperator(n,n;kmin=6.,kmax=8.,allow_zonal=false)
m=build_model(o;f=1.);u,v=random_seed(o,seed);seed!(m,u,v)
b=Audit(o;dt)
config=Dict("grid"=>[n,n],"domain"=>[2pi,2pi],"geometry"=>"periodic horizontal 2D",
    "seed"=>seed,"dt"=>dt,"target_growth_times"=>target,"target_iterations"=>steps,
    "f"=>1.0,"planetary_beta"=>0.0,"amplitude"=>b.amplitude,"viscosity"=>b.nu,
    "cubic_velocity_drag"=>b.alpha,"forcing_band"=>[6.,8.],"allow_zonal_forcing"=>false,
    "remove_mean_drag"=>true,
    "sigma_max_injection"=>sigma,"initial_rms_velocity"=>0.02,
    "physical_time_mapping"=>"none; dimensionless time, not Jupiter days",
    "oceananigans_version"=>string(pkgversion(Oceananigans)),"compute_location"=>"Unity",
    "wall_time_limit_seconds"=>2400,"scientific_scope"=>"Exploratory instability-forced 2D pilot; not an exact van Kan reproduction",
    "source_sha256"=>Dict(basename(p)=>bytes2hex(sha256(read(p))) for p in
        [@__FILE__,joinpath(@__DIR__,"audit.jl"),joinpath(@__DIR__,"../instability2d_v1/core.jl")]))
open(joinpath(out,"config.toml"),"w") do io;TOML.print(io,config;sorted=true);end
open(joinpath(out,"statistics.csv"),"w") do io
    println(io,"iteration,t,growth_times,energy,enstrophy,skewness,zeta_min,zeta_max,low_k_energy_fraction,advection_M3,instability_M3,viscous_M3,drag_M3,numerical_advection_fraction")
end
open(joinpath(out,"frames.csv"),"w") do io;println(io,"iteration,t,growth_times,file");end
started=time();lastframe=Ref(-1);laststat=Ref(-1)
function output!(sim;force=false)
    m=sim.model;finish!(b,m);it=m.clock.iteration
    if (it%50==0 || force) && laststat[]!=it
        u,v=state(m);q=vorticity(o,u,v);qr=sqrt(mean(q.^2));energy=mean(u.^2 .+ v.^2)/2
        qh=FFTW.fft(q);invlam=map(x->x==0 ? 0.0 : 1/x,o.lambda)
        ek=abs2.(qh).*invlam;low=sum(ek[o.k2.<=4])/max(sum(ek),1e-30)
        ratio=b.absolute_third[1]/max(sum(b.absolute_third[2:4]),1e-30)
        open(joinpath(out,"statistics.csv"),"a") do io
            println(io,join((it,m.clock.time,sigma*m.clock.time,energy,qr^2/2,mean(q.^3)/max(qr^3,1e-30),
                minimum(q),maximum(q),low,b.third...,ratio),','))
        end
        laststat[]=it
    end
    if (it%200==0 || force) && lastframe[]!=it
        u,v=state(m);q=vorticity(o,u,v);name=@sprintf("zeta_%07d.f32",it)
        open(joinpath(out,name),"w") do io;write(io,Float32.(q));end
        open(joinpath(out,"frames.csv"),"a") do io;println(io,join((it,m.clock.time,sigma*m.clock.time,name),','));end
        lastframe[]=it
        println(@sprintf("PROGRESS iteration=%d growth_times=%.2f elapsed=%.1fs CFL=%.3f budget_error=%.3g",it,sigma*m.clock.time,time()-started,b.max_cfl,b.max_increment_error));flush(stdout)
    end
    if it>0 && it%2000==0
        u,v=state(m);jldsave(joinpath(out,"checkpoint.tmp.jld2");u,v,t=Float64(m.clock.time),iteration=it,config)
        mv(joinpath(out,"checkpoint.tmp.jld2"),joinpath(out,"checkpoint.jld2");force=true)
    end
end
sim=Simulation(m;Δt=dt,stop_iteration=steps,wall_time_limit=2400.,verbose=false)
attach_audit!(sim,b)
add_callback!(sim,output!,IterationInterval(1);name=:pilot_output)
run!(sim);output!(sim;force=true)
u,v=state(m);jldsave(joinpath(out,"checkpoint.jld2");u,v,t=Float64(m.clock.time),iteration=m.clock.iteration,config)
result=report(b);result["status"]=m.clock.iteration>=steps ? "COMPLETED" : "PARTIAL_WALL_TIME_LIMIT"
result["growth_times_completed"]=sigma*m.clock.time;result["elapsed_seconds"]=time()-started
result["scope"]="Exploratory pilot; no sign preference or inverse-cascade claim from a single realization"
open(joinpath(out,"review.toml"),"w") do io;TOML.print(io,result;sorted=true);end
println("PILOT_FINISHED");TOML.print(stdout,result;sorted=true)

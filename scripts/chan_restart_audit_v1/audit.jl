include(joinpath(@__DIR__, "../chan_clean_anelastic6000_v1/model.jl"))
function main()
    out=ARGS[1];@assert !ispath(out);mkpath(out)
    for name in ("chan_clean_anelastic6000_v1", "chan_theta_diagnostic4000_v1")
        cfg=TOML.parsefile(joinpath(@__DIR__,"..",name,"config.toml"))
        h=TOML.parsefile(joinpath("runs",replace(name,"_"=>"-"),"HEALTH_PASS.toml"))
        m,p=build_model(cfg;arch=CPU())
        sim=Simulation(m;Δt=cfg["initial_dt"],stop_time=cfg["restart_time"])
        set!(sim;checkpoint=h["restart_checkpoint"])
        @assert abs(m.clock.time-cfg["restart_time"])<1e-9
        Oceananigans.TimeSteppers.update_state!(m)
        a=arrays(m);dv=p.dx*p.dy*p.dz
        rx=.5 .* (a.rho.+circshift(a.rho,(1,0,0)));ry=.5 .* (a.rho.+circshift(a.rho,(0,1,0)))
        rz=.5 .* (a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
        ke=.5dv*(sum(rx.*a.uf.^2)+sum(ry.*a.vf.^2)+sum(rz.*a.wf[:,:,2:end-1].^2))
        z=[zcell(k,p) for k=1:p.nz]
        e=dv*sum(a.rho.*(CV.*a.T.+GRAV.*reshape(z,1,1,:)))+ke
        checks=Dict("kinetic_energy"=>ke,"energy_proxy"=>e,"min_density"=>minimum(a.rho),"min_temperature"=>minimum(a.T))
        errors=Dict(k=>abs(v-h["donor_"*k])/max(abs(h["donor_"*k]),1.) for (k,v) in checks)
        err=maximum(values(errors));@assert err<1e-10
        result=Dict("status"=>"PASS","time"=>Float64(m.clock.time),"iteration"=>m.clock.iteration,"max_relative_error"=>err,"values"=>checks,"relative_errors"=>errors,"note"=>"Independent CPU checkpoint load; derived thermodynamics refreshed; no integration performed")
        open(joinpath(out,name*".toml"),"w") do io;TOML.print(io,result);end
        println(name," NATIVE_CHECKPOINT_AUDIT_PASS ",err);flush(stdout)
        m=nothing;sim=nothing;a=nothing;GC.gc()
    end
end
main()

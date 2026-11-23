function save_bad(m,out,reason)
    path=joinpath(out,"FIRST_BAD_STATE.jld2")
    isfile(path) && return
    jldsave(path;t=Float64(m.clock.time),iteration=m.clock.iteration,rho=Array(interior(m.dynamics.dry_density)),rhotheta=Array(interior(m.formulation.potential_temperature_density)),T=Array(interior(m.temperature)),P=Array(interior(m.dynamics.pressure)),u=Array(interior(m.velocities.u)),v=Array(interior(m.velocities.v)),w=Array(interior(m.velocities.w)),nu=Array(interior(m.closure_fields.νₑ)))
    write_toml(joinpath(out,"FIRST_BAD_STATE.toml"),Dict("time"=>Float64(m.clock.time),"iteration"=>m.clock.iteration,"reason"=>reason))
end

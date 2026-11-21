function reflect_x!(m,cfg,out)
    nx=cfg["nx"]
    report=Dict{String,Any}("operation"=>"x -> -x; u -> -u; v,w and scalars reflected without sign change; f reversed separately")
    for (key,field) in pairs(Oceananigans.prognostic_fields(m))
        a=Array(interior(field));name=string(key)
        q=name=="ρu" ? -a[mod1.(2 .- (1:nx),nx),:,:] : reverse(a;dims=1)
        set!(field,q)
        back=name=="ρu" ? -q[mod1.(2 .- (1:nx),nx),:,:] : reverse(q;dims=1)
        back==a || error("Reflection does not invert exactly: $name")
    end
    set!(m.dynamics.pressure_anomaly,reverse(Array(interior(m.dynamics.pressure_anomaly));dims=1))
    report["status"]="PASS"
    open(joinpath(out,"reflection.toml"),"w") do io;TOML.print(io,report);end
end

# All manipulations act on native staggered momentum, never interpolated movie fields.
function reverse_barotropic!(m,rho,cfg,out)
    fields=Oceananigans.prognostic_fields(m)
    old=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(fields))
    h=reshape(rho,1,1,:);weights=h./sum(rho)
    reports=Dict{String,Any}("operation"=>"Reverse density-weighted depth-mean horizontal velocity; preserve 3D deviations, w and thermodynamics")
    for (key,field) in pairs(fields)
        name=string(key)
        if name in ["ρu","ρv"]
            velocity=old[name]./h;bar=sum(weights.*velocity;dims=3)
            q=old[name].-2 .*h.*bar
            newbar=sum(weights.*(q./h);dims=3)
            err=maximum(abs,newbar.+bar)
            err<1e-11 || error("Barotropic reversal failed")
            reports[name*"_barotropic_reversal_error_m_s"]=err
            reports[name*"_eddy_preservation_error_m_s"]=maximum(abs,(q./h.-newbar).-(velocity.-bar))
            set!(field,q)
        end
    end
    new=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(fields))
    horizontal_ke(a)=sum(a["ρu"].^2 ./h)+sum(a["ρv"].^2 ./h)
    reports["horizontal_kinetic_energy_relative_change"]=horizontal_ke(new)/horizontal_ke(old)-1
    abs(reports["horizontal_kinetic_energy_relative_change"])<1e-12 || error("Reversal changed horizontal kinetic energy")
    dx,dy,dz=cfg["lx"]/cfg["nx"],cfg["ly"]/cfg["ny"],cfg["height"]/cfg["nz"]
    div(a)=(circshift(a["ρu"],(-1,0,0)).-a["ρu"])./dx.+diff(a["ρv"];dims=2)./dy.+diff(a["ρw"];dims=3)./dz
    reports["mass_divergence_before"]=maximum(abs,div(old));reports["mass_divergence_after"]=maximum(abs,div(new))
    reports["mass_divergence_after"]<1e-11 || error("Reversed state is not divergence free")
    for key in keys(old)
        key in ["ρu","ρv"] && continue
        old[key]==new[key] || error("Reversal modified $key")
    end
    reports["pressure_note"]="Original pressure is retained as an initial guess; the anelastic projection adjusts it on the first step. Compare with an unchanged restart."
    reports["status"]="PASS"
    open(joinpath(out,"reversal_validation.toml"),"w") do io;TOML.print(io,reports);end
end

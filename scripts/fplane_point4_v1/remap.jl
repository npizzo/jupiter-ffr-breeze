# Conservative horizontal refinement of staggered mass fluxes on Unity.
# D_f P_F = P_C D_c, so the same linear center prolongation in transverse
# directions preserves discrete anelastic continuity without a new projection.
function prolong_center(a,dim;periodic=false,limited=false)
    n=size(a,dim);sz=collect(size(a));sz[dim]=2n;b=Array{eltype(a)}(undef,sz...)
    slab(j)=selectdim(a,dim,periodic ? mod1(j,n) : clamp(j,1,n))
    for j=1:n
        q=slab(j);qm=slab(j-1);qp=slab(j+1)
        delta = if limited
            left=q.-qm;right=qp.-q
            ifelse.(left.*right.>0,sign.(left).*min.(abs.(left),abs.(right))./4,zero(eltype(a)))
        else
            (qp.-qm)./8
        end
        selectdim(b,dim,2j-1).=q.-delta
        selectdim(b,dim,2j).=q.+delta
    end
    b
end
function prolong_faces(a,dim;periodic=false)
    n=size(a,dim)-(periodic ? 0 : 1);sz=collect(size(a));sz[dim]=2n+(periodic ? 0 : 1)
    b=Array{eltype(a)}(undef,sz...)
    function slab(j)
        periodic && return selectdim(a,dim,mod1(j,n))
        j==0 && return -selectdim(a,dim,2)
        j==n+2 && return -selectdim(a,dim,n)
        selectdim(a,dim,j)
    end
    for j=1:n
        selectdim(b,dim,2j-1).=slab(j)
        selectdim(b,dim,2j).=(slab(j).+slab(j+1))./2 .- (slab(j+2).-slab(j+1).-slab(j).+slab(j-1))./16
    end
    periodic || (selectdim(b,dim,2n+1).=slab(n+1))
    b
end
function average_pairs(a,dim)
    n=size(a,dim);iseven(n) || error("Pair average requires even size")
    (selectdim(a,dim,1:2:n).+selectdim(a,dim,2:2:n))./2
end
function refine_field(a,key)
    if key=="ρu"
        prolong_center(prolong_faces(a,1;periodic=true),2)
    elseif key=="ρv"
        prolong_center(prolong_faces(a,2),1;periodic=true)
    else
        limited=key!="ρw" && key!="pressure"
        prolong_center(prolong_center(a,1;periodic=true,limited),2;limited)
    end
end
function restrict_field(a,key)
    if key=="ρu"
        average_pairs(a[1:2:end,:,:],2)
    elseif key=="ρv"
        average_pairs(a[:,1:2:end,:],1)
    else
        average_pairs(average_pairs(a,1),2)
    end
end
mass_divergence(M,dx,dy,dz)=(circshift(M["ρu"],(-1,0,0)).-M["ρu"])./dx .+ diff(M["ρv"];dims=2)./dy .+ diff(M["ρw"];dims=3)./dz
function export_bundle(m,sim,controller,path,cfg,pickup)
    isfile(path) && error("Refusing to overwrite restart bundle")
    prognostic=Dict(string(k)=>Array(interior(v)) for (k,v) in pairs(Oceananigans.prognostic_fields(m)))
    all(k->haskey(prognostic,k),["ρu","ρv","ρw","ρθ","dye"]) || error("Unexpected model fields")
    jldsave(path;prognostic,pressure=Array(interior(m.dynamics.pressure_anomaly)),
        clock_state=Oceananigans.prognostic_state(m.clock),controller_state=Oceananigans.prognostic_state(controller),
        dt=sim.Δt,density=Array(interior(AM.continuity_density(m.dynamics))),
        grid_size=[cfg["nx"],cfg["ny"],cfg["nz"]],domain=[cfg["lx"],cfg["ly"],cfg["height"]],
        source_checkpoint=pickup,source_checkpoint_sha256=bytes2hex(open(sha256,pickup)))
    println("EXPORTED_DAY90_BUNDLE ",path);flush(stdout)
end
function restore_bundle!(m,sim,controller,saved_dt,cfg,out)
    b=load(cfg["restart_bundle"]);oldsize=b["grid_size"];newsize=[cfg["nx"],cfg["ny"],cfg["nz"]]
    oldsize==[256,256,64] || error("Unexpected source resolution")
    newsize in ([256,256,64],[512,512,64]) || error("Only same-grid or horizontal factor-two refinement")
    b["domain"]==[cfg["lx"],cfg["ly"],cfg["height"]] || error("Domain changed")
    rho=Array(interior(AM.continuity_density(m.dynamics)))
    rho==b["density"] || error("Reference density changed")
    fine=newsize[1]==2oldsize[1]
    old=b["prognostic"];reports=Dict{String,Any}();transferred=Dict{String,Array{Float64,3}}()
    for (key,field) in pairs(Oceananigans.prognostic_fields(m))
        name=string(key);a=old[name];q=fine ? refine_field(a,name) : a
        size(q)==size(interior(field)) || error("Wrong staggered shape for $name")
        all(isfinite,q) || error("Nonfinite remap")
        back=fine ? restrict_field(q,name) : q
        err=maximum(abs,back.-a);scale=max(maximum(abs,a),1.0)
        err<=1e-12scale || error("Conservative restriction failed for $name: $err")
        if name in ["ρθ","ρqᵛ","dye"]
            minimum(q)>=minimum(a)-1e-12scale && maximum(q)<=maximum(a)+1e-12scale || error("Scalar remap creates extrema")
        end
        reports[name*"_restriction_error"]=err
        reports[name*"_relative_l2_change_on_restriction"]=sqrt(sum(abs2,back.-a)/max(sum(abs2,a),1e-30))
        set!(field,q)
        name in ["ρu","ρv","ρw"] && (transferred[name]=q)
    end
    dx,dy,dz=b["domain"]./oldsize
    olddiv=mass_divergence(old,dx,dy,dz)
    newdiv=mass_divergence(transferred,cfg["lx"]/cfg["nx"],cfg["ly"]/cfg["ny"],dz)
    expected=fine ? prolong_center(prolong_center(olddiv,1;periodic=true),2) : olddiv
    commutation=maximum(abs,newdiv.-expected)
    commutation<1e-12 || error("Divergence-preserving remap failed: $commutation")
    maximum(abs,newdiv)<1e-11 || error("Restart is not discretely divergence-free")
    maximum(abs,transferred["ρv"][:,[1,end],:])==0 || error("Meridional boundary flux changed")
    maximum(abs,transferred["ρw"][:,:,[1,end]])==0 || error("Vertical boundary flux changed")
    set!(m.dynamics.pressure_anomaly,fine ? refine_field(b["pressure"],"pressure") : b["pressure"])
    Oceananigans.restore_prognostic_state!(m.clock,b["clock_state"])
    Oceananigans.restore_prognostic_state!(controller,b["controller_state"])
    saved_dt.value=min(b["dt"],cfg["max_dt"]);sim.Δt=saved_dt.value
    m.clock.time==90*86400.0 && m.clock.stage==1 || error("Wrong restart time or stage")
    merge!(reports,Dict("status"=>"PASS","source_checkpoint_sha256"=>b["source_checkpoint_sha256"],
        "source_grid"=>oldsize,"target_grid"=>newsize,"mass_divergence_commutation_error"=>commutation,
        "source_max_mass_divergence"=>maximum(abs,olddiv),"target_max_mass_divergence"=>maximum(abs,newdiv),
        "momentum_method"=>"Linear conservative center prolongation paired with commuting normal-face prolongation",
        "scalar_method"=>"Conservative minmod-limited horizontal slopes","pressure_projection_applied"=>false,
        "timestep_seconds"=>sim.Δt,"compute_location"=>"Unity"))
    open(joinpath(out,"remap_validation.toml"),"w") do io;TOML.print(io,reports);end
    nothing
end

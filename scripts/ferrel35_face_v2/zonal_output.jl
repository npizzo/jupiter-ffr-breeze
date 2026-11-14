# Read-only centered covariance diagnostics and native meridional mass fluxes.
function zonal_output(m, u, v, w, theta, dye, dy, dz)
    zm(a)=dropdims(mean(a;dims=1);dims=1)
    zonal_u,zonal_v,zonal_w=zm(u),zm(v),zm(w)
    zonal_theta,zonal_tracer=zm(theta),zm(dye)
    eddy_uv=zm(u.*v).-zonal_u.*zonal_v
    eddy_uw=zm(u.*w).-zonal_u.*zonal_w
    eddy_vtheta=zm(v.*theta).-zonal_v.*zonal_theta
    eddy_wtheta=zm(w.*theta).-zonal_w.*zonal_theta
    eddy_vtracer=zm(v.*dye).-zonal_v.*zonal_tracer
    eddy_wtracer=zm(w.*dye).-zonal_w.*zonal_tracer
    # Read the actual native mass fluxes; do not approximate face density.
    native_rhov=zm(Array(interior(m.momentum.ρv)))
    native_rhow=zm(Array(interior(m.momentum.ρw)))
    ny,nz=size(zonal_u)
    size(native_rhov)==(ny+1,nz) && size(native_rhow)==(ny,nz+1) || error("Native flux staggering")
    mass_streamfunction=zeros(ny+1,nz+1)
    mass_streamfunction[:,2:end].=-cumsum(native_rhov;dims=2).*dz
    reconstructed_w=diff(mass_streamfunction;dims=1)./dy
    mass_div=diff(native_rhov;dims=1)./dy.+diff(native_rhow;dims=2)./dz
    return (;zonal_u,zonal_v,zonal_w,zonal_theta,zonal_tracer,
        eddy_uv,eddy_uw,eddy_vtheta,eddy_wtheta,eddy_vtracer,eddy_wtracer,
        native_rhov,native_rhow,mass_streamfunction,
        mass_divergence_max=maximum(abs,mass_div),
        mass_flux_reconstruction_max=maximum(abs,reconstructed_w.-native_rhow),
        top_streamfunction_max=maximum(abs,mass_streamfunction[:,end]),
        zonal_diagnostic_note="Zonal means and centered resolved covariances; not an exact discrete momentum budget. Native mass streamfunction: rho*v=-dPsi/dz, rho*w=dPsi/dy, bottom Psi=0, units kg/(m*s).")
end

function check_zonal_frames(out)
    m=TOML.parsefile(joinpath(out,"metadata.toml"));ny,nz=m["ny"],m["nz"]
    for p in filter(p->startswith(basename(p),"frame_") && endswith(p,".jld2"),readdir(out;join=true))
        jldopen(p) do d
            for k in ("zonal_u","zonal_v","zonal_w","zonal_theta","zonal_tracer","eddy_uv","eddy_uw","eddy_vtheta","eddy_wtheta","eddy_vtracer","eddy_wtracer")
                @assert size(d[k])==(ny,nz) && all(isfinite,d[k])
            end
            @assert size(d["native_rhov"])==(ny+1,nz)
            @assert size(d["native_rhow"])==(ny,nz+1)
            @assert size(d["mass_streamfunction"])==(ny+1,nz+1)
            for k in ("native_rhov","native_rhow","mass_streamfunction")
                @assert all(isfinite,d[k])
            end
        end
    end
end

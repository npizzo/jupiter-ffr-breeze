"""A wall-compatible barotropic cosine jet with discrete beta-geostrophic pressure."""
function balanced_cosine_jet(cfg,f)
    ny=cfg["ny"];ly=cfg["ly"];dy=ly/ny
    wavelength=cfg["jet_wavelength"];speed=cfg["initial_jet_speed"]
    cycles=ly/wavelength
    isapprox(cycles,round(cycles);atol=1e-12) || error("Cosine jet must have an integer number of periods in channel")
    yc=((1:ny).-0.5).*dy.-ly/2
    u=speed.*cos.(2pi.*(yc.+ly/2)./wavelength)
    yv=(1:ny-1).*dy.-ly/2
    force=(f.+cfg["beta"].*yv).*(u[1:end-1].+u[2:end])./2
    pressure=vcat(0.0,-cumsum(force).*dy)
    pressure.-=mean(pressure)
    residual=maximum(abs,(pressure[2:end].-pressure[1:end-1])./dy.+force)
    residual<1e-12*max(abs(f*speed),1e-15) || error("Discrete jet pressure balance failed")
    ratio=abs(speed)*(2pi/wavelength)^2/max(abs(cfg["beta"]),eps())
    return u,pressure,Dict{String,Any}("initial_jet_model"=>"Barotropic cosine, no maintained force", 
        "initial_jet_pressure_residual"=>residual,"initial_jet_curvature_over_beta"=>ratio,
        "initial_jet_wavelength_m"=>wavelength,"initial_jet_speed_m_s"=>speed,
        "initial_jet_pressure_balance"=>"Discrete integral of -f(y_face) times centered U at v faces")
end

# Centered products are for interpretation, not the exact WENO transport.
# One-sided derivatives at y/z boundaries are explicit in the metadata.
function centered_derivative(a,h,dim)
    out=similar(a)
    if dim==1
        out[2:end-1,:].=(a[3:end,:].-a[1:end-2,:])./(2h)
        out[1,:].=(a[2,:].-a[1,:])./h;out[end,:].=(a[end,:].-a[end-1,:])./h
    else
        out[:,2:end-1].=(a[:,3:end].-a[:,1:end-2])./(2h)
        out[:,1].=(a[:,2].-a[:,1])./h;out[:,end].=(a[:,end].-a[:,end-1])./h
    end
    return out
end

function centered_flux_interpretation(u,v,w,uv,uw,rho,dy,dz,native_advection)
    density=reshape(rho,1,:)
    mean_horizontal=-centered_derivative(density.*u.*v,dy,1)
    mean_vertical=-centered_derivative(density.*u.*w,dz,2)
    eddy_horizontal=-centered_derivative(uv,dy,1)
    eddy_vertical=-centered_derivative(uw,dz,2)
    remainder=native_advection.-mean_horizontal.-mean_vertical.-eddy_horizontal.-eddy_vertical
    # Positive means eddy covariance stresses feed the zonal mean kinetic energy.
    transfer_y=uv.*centered_derivative(u,dy,1)
    transfer_z=uw.*centered_derivative(u,dz,2)
    return (;centered_mean_advection_y=mean_horizontal,centered_mean_advection_z=mean_vertical,
        centered_eddy_advection_y=eddy_horizontal,centered_eddy_advection_z=eddy_vertical,
        centered_native_advection_remainder=remainder,
        eddy_to_mean_ke_transfer_y=transfer_y,eddy_to_mean_ke_transfer_z=transfer_z,
        eddy_to_mean_ke_transfer_W_m2=dz*sum(transfer_y.+transfer_z)/size(u,1),
        centered_eddy_force_work_W_m2=dz*sum(u.*(eddy_horizontal.+eddy_vertical))/size(u,1),
        mean_ke_J_m2=dz*sum(density.*u.^2)/2size(u,1),
        centered_flux_definition="Interpolated covariances; centered interior/one-sided boundary derivatives; remainder against native WENO advection retained. Positive covariance*mean shear feeds mean KE.")
end

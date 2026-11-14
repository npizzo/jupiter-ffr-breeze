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

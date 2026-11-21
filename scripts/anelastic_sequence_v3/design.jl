include("model.jl")
root=ARGS[1];mkpath(root)
c=TOML.parsefile("configs/anelastic-sequence-v3/baseline.toml")
reports=Dict{String,Any}()
for (label,caplog) in (("weak",log(1.01)),("strong",.25))
    cc=merge(c,Dict("nx"=>8,"ny"=>8,"cap_log_theta"=>caplog))
    m,p=build(cc);a=fields(m);z=[(k-.5)*p.dz for k=1:p.nz]
    rho=vec(a.rho);thref=theta_reference.(z,Ref(p));th=vec(mean(a.theta;dims=(1,2)))
    fluxref=-p.kappa.*(0.5 .*(rho[1:end-1].+rho[2:end])).*diff(thref)./p.dz
    N2=p.g.*diff(log.(thref))./p.dz
    capcool=sum(rho.*shape.(z,Ref(p)).*(th.-thref))/p.tau*p.dz
    B=p.g*p.flux;U=(B*p.dc)^(1/3);E=p.nu/(abs(p.f)*p.dc^2)
    report=Dict("cap_log_theta"=>caplog,"theta_top_over_base"=>exp(caplog),"Tmin_reference"=>minimum(Array(interior(p.Tref))),"density_contrast_convective"=>rho[1]/rho[argmin(abs.(z.-p.dc))],"density_contrast_full"=>rho[1]/rho[end],"N_over_f_max"=>sqrt(maximum(N2))/abs(p.f),"N_over_f_rms_cap"=>sqrt(p.g*caplog/(p.height-p.dc))/abs(p.f),"cap_min_tau_f"=>p.tau*abs(p.f),"Pr"=>p.nu/p.kappa,"Ek"=>E,"flux_velocity_scale_bottom"=>U,"flux_Ro_depth_estimate"=>U/(abs(p.f)*p.dc),"flux_Re_depth_estimate"=>U*p.dc/p.nu,"onset_wavelength_screen"=>4.8*E^(1/3)*p.dc,"onset_wavelength_cells_screen"=>4.8*E^(1/3)*p.dc/(p.length/c["nx"]),"bottom_layer_cells"=>p.bl/p.dz,"bottom_layer_theta_excess"=>p.flux*p.bl/p.kappa,"bottom_layer_diffusive_time"=>p.bl^2/p.kappa,"initial_cap_cooling_over_bottom_input"=>capcool/p.flux,"reference_downward_diffusion_peak_over_input"=>max(0.,-minimum(fluxref))/p.flux,"note"=>"Reference-state diffusion may dominate cap thermal balance. This is a source calculation, not developed convection or a launch approval.")
    reports[label]=report
    println("DESIGN ",label," N/f=",report["N_over_f_rms_cap"]," reference_downflux/F=",report["reference_downward_diffusion_peak_over_input"]," lambda_cells=",report["onset_wavelength_cells_screen"]);flush(stdout)
end
record(joinpath(root,"DESIGN.toml"),reports)

using Random, Statistics

"""Full-depth random streamfunction, discretely nondivergent and wall-compatible.

u = +ψ_y, v = -ψ_x. The cell pressure estimate is -f0 ψ, so it balances
the linear f0 Coriolis term with centered C-grid interpolation. This is not
an exact nonlinear equilibrium or an exact beta-plane geostrophic balance.
"""
function initial_large_eddies(cfg, rho, f)
    nx,ny,nz=cfg["nx"],cfg["ny"],cfg["nz"]
    lx,ly,h=cfg["lx"],cfg["ly"],cfg["height"]
    dx,dy,dz=lx/nx,ly/ny,h/nz
    wavelength=cfg["turbulence_wavelength"]; k0=2pi/wavelength
    cfg["turbulence_vertical_ratio"]==0 || error("Barotropic seed requires w=0")
    wavelength/max(dx,dy)/1.25 >= 16 || error("Shortest seeded wavelength is under-resolved")
    length(rho)==nz && all(rho .> 0) || error("Invalid reference density")
    target=cfg["turbulence_rossby"]*abs(f)/k0
    target>0 || error("Nonzero Rossby number and rotation required")
    rng=MersenneTwister(cfg["turbulence_seed"])
    # Corners: periodic x faces and both y walls; omit zonal-mean modes.
    psi=zeros(nx,ny+1); modes=Tuple{Int,Int,Float64}[]
    for m in 1:ceil(Int,1.25lx/wavelength), n in 1:ceil(Int,2.5ly/wavelength)
        kh=hypot(2pi*m/lx,pi*n/ly)
        0.75k0 <= kh <= 1.25k0 || continue
        a=randn(rng); phase=2pi*rand(rng)
        psi .+= a .* reshape(cos.(2pi*m.*(0:nx-1)./nx .+ phase),nx,1) .* reshape(sinpi.(n.*(0:ny)./ny),1,ny+1)
        push!(modes,(m,n,kh))
    end
    isempty(modes) && error("No random modes in requested wavelength band")
    psi[:,[1,end]].=0
    u=(psi[:,2:end].-psi[:,1:end-1])./dy
    v=-(circshift(psi,(-1,0)).-psi)./dx
    rms=sqrt(mean(u.^2 .+ (v[:,1:end-1].^2 .+v[:,2:end].^2)./2))
    scale=target/rms; psi .*=scale; u .*=scale; v .*=scale
    div=(circshift(u,(-1,0)).-u)./dx .+ (v[:,2:end].-v[:,1:end-1])./dy
    relative_div=maximum(abs,div)/(target/min(dx,dy))
    relative_div<1e-12 || error("Seed violates discrete divergence constraint")
    maximum(abs,v[:,[1,end]])==0 || error("Nonzero wall-normal flow")
    pc=-f .* (psi[:,1:end-1].+psi[:,2:end].+circshift(psi[:,1:end-1],(-1,0)).+circshift(psi[:,2:end],(-1,0)))./4
    # Check the linear f0 balance at staggered horizontal-velocity points.
    v_at_u=(v[:,1:end-1].+v[:,2:end].+circshift(v[:,1:end-1],(1,0)).+circshift(v[:,2:end],(1,0)))./4
    px=(pc.-circshift(pc,(1,0)))./dx
    u_at_v=(u[:,1:end-1].+u[:,2:end].+circshift(u[:,1:end-1],(-1,0)).+circshift(u[:,2:end],(-1,0)))./4
    py=(pc[:,2:end].-pc[:,1:end-1])./dy
    residual=max(maximum(abs,px.-f.*v_at_u),maximum(abs,py.+f.*u_at_v))/(abs(f)*target)
    residual<1e-12 || error("f0-geostrophic seed pressure is inconsistent")
    info=Dict{String,Any}(
        "initial_velocity_model"=>"Full-depth barotropic random streamfunction; no zonal-mean random mode; zero w",
        "initial_horizontal_rms_lower_m_s"=>target,"initial_vertical_rms_lower_m_s"=>0.0,
        "initial_eddy_rossby"=>cfg["turbulence_rossby"],"initial_eddy_length_m"=>1/k0,
        "initial_mass_divergence_relative"=>relative_div,"initial_horizontal_mode_count"=>length(modes),
        "initial_seed_wavelength_min_m"=>minimum(2pi/t[3] for t in modes),
        "initial_seed_wavelength_max_m"=>maximum(2pi/t[3] for t in modes),
        "initial_random_geostrophic_residual_f0"=>residual,
        "initial_pressure_balance"=>"Centered f0-geostrophic pressure estimate; beta and nonlinear adjustment retained",
        "initial_beta_fractional_f_variation"=>abs(cfg["beta"])*ly/(2abs(f)),
        "initial_random_vertical_structure"=>"Identical horizontal velocity at every height, including stable cap")
    return repeat(u,1,1,nz),repeat(v,1,1,nz),zeros(nx,ny,nz+1),repeat(pc,1,1,nz),info
end

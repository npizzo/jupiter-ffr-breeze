using Random, Statistics

"""Band-limited toroidal/poloidal mass fluxes on a uniform staggered grid."""
function initial_channel_turbulence(cfg, rho, f)
    nx, ny, nz = cfg["nx"], 2cfg["ny"], cfg["nz"]
    half = ny÷2
    lx, ly, h, depth = cfg["lx"], 2cfg["ly"], cfg["height"], cfg["cap_height"]
    dx, dy, dz = lx/nx, ly/ny, h/nz
    wavelength = cfg["turbulence_wavelength"]
    length_scale = wavelength/(2pi)
    target = cfg["turbulence_rossby"]*abs(f)*length_scale
    target > 0 || error("Turbulence normalization requires nonzero rotation and target Rossby number")
    wavelength/min(dx,dy) >= 16 || error("Initial eddies require at least 16 cells per wavelength")
    rng = MersenneTwister(cfg["turbulence_seed"])
    zc = ((1:nz).-0.5).*dz; zf = (0:nz).*dz
    # Psi lives on vertical edges; Phi lives on vertical-velocity faces.
    psi = zeros(nx,ny,nz); phi = zeros(nx,ny,nz+1)
    mode_count = 0
    k0 = 2pi/wavelength
    for m in 0:ceil(Int,1.25lx/wavelength), n in 0:ceil(Int,1.25ly/wavelength)
        (m > 0 || n > 0) || continue
        kh = hypot(2pi*m/lx,2pi*n/ly)
        0.75k0 <= kh <= 1.25k0 || continue
        mode_count += 1
        for vertical_mode in 0:2
            a, b = randn(rng), randn(rng)
            phase_a, phase_b = 2pi*rand(rng), 2pi*rand(rng)
            horizontal_a = [cos(2pi*m*(i-1)/nx+phase_a)*sin(2pi*n*(j-1)/ny) for i=1:nx,j=1:ny]
            horizontal_b = [cos(2pi*m*(i-0.5)/nx+phase_b)*cos(2pi*n*(j-0.5)/ny) for i=1:nx,j=1:ny]
            shape_a = [0 < z < depth ? sinpi(z/depth)^2*cospi(vertical_mode*z/depth) : 0.0 for z in zc]
            shape_b = [0 < z < depth ? sinpi(z/depth)^3*cospi(vertical_mode*z/depth) : 0.0 for z in zf]
            psi .+= a.*reshape(horizontal_a,nx,ny,1).*reshape(shape_a,1,1,nz)
            phi .+= b.*reshape(horizontal_b,nx,ny,1).*reshape(shape_b,1,1,nz+1)
        end
    end
    mode_count > 0 || error("No resolved Fourier modes in turbulence band")
    # Exact odd/even extension gives no-normal-flow and free-slip walls.
    psi[:,1,:] .= 0; psi[:,half+1,:] .= 0
    psi[:,half+2:end,:] .= -psi[:,half:-1:2,:]
    phi[:,half+1:end,:] .= phi[:,half:-1:1,:]
    rho_c = reshape(rho,1,1,nz)
    rho_f = reshape(vcat(rho[1],(rho[1:end-1].+rho[2:end])./2,rho[end]),1,1,nz+1)
    ut = (circshift(psi,(0,-1,0)).-psi)./dy./rho_c
    vt = -(circshift(psi,(-1,0,0)).-psi)./dx./rho_c
    gx = (phi.-circshift(phi,(1,0,0)))./dx
    gy = (phi.-circshift(phi,(0,1,0)))./dy
    up = (gx[:,:,2:end].-gx[:,:,1:end-1])./dz./rho_c
    vp = (gy[:,:,2:end].-gy[:,:,1:end-1])./dz./rho_c
    wp = -((circshift(gx,(-1,0,0)).-gx)./dx.+(circshift(gy,(0,-1,0)).-gy)./dy)./rho_f
    lower = zc .< depth
    wc = (wp[:,:,1:end-1].+wp[:,:,2:end])./2
    vertical_target = cfg["turbulence_vertical_ratio"]*target
    ap = vertical_target/sqrt(mean(abs2,wc[:,:,lower]))
    # Solve the quadratic including the small discrete toroidal/poloidal cross term.
    aa = mean(ut[:,:,lower].^2 .+vt[:,:,lower].^2)
    bb = 2ap*mean(ut[:,:,lower].*up[:,:,lower].+vt[:,:,lower].*vp[:,:,lower])
    cc = ap^2*mean(up[:,:,lower].^2 .+vp[:,:,lower].^2)-target^2
    cc < 0 || error("Requested vertical component is too large for the horizontal RMS target")
    at = (-bb+sqrt(bb^2-4aa*cc))/(2aa)
    u = at.*ut.+ap.*up; v = at.*vt.+ap.*vp; w = ap.*wp
    mass_u, mass_v, mass_w = rho_c.*u, rho_c.*v, rho_f.*w
    div = (circshift(mass_u,(-1,0,0)).-mass_u)./dx.+(circshift(mass_v,(0,-1,0)).-mass_v)./dy.+(mass_w[:,:,2:end].-mass_w[:,:,1:end-1])./dz
    relative_div = maximum(abs,div)/(maximum(rho)*target/min(dx,dy,dz))
    relative_div < 1e-10 || error("Initial mass-flux divergence exceeds tolerance")
    maximum(abs,w[:,:,[1,end]]) == 0 || error("Nonzero vertical boundary velocity")
    measured = sqrt(mean(u[:,:,lower].^2 .+v[:,:,lower].^2))
    measured_w = sqrt(mean(abs2,(w[:,:,1:end-1].+w[:,:,2:end])[:,:,lower]./2))
    info = Dict("initial_velocity_model"=>"band-limited toroidal and poloidal mass flux",
        "initial_eddy_length_m"=>length_scale,"initial_horizontal_rms_lower_m_s"=>measured,
        "initial_vertical_rms_lower_m_s"=>measured_w,"initial_eddy_rossby"=>measured/(abs(f)*length_scale),
        "initial_mass_divergence_relative"=>relative_div,"initial_horizontal_mode_count"=>mode_count)
    maximum(abs,v[:,[1,half+1],:]) == 0 || error("Nonzero meridional wall velocity")
    info["initial_velocity_model"] = "wall-compatible toroidal/poloidal mass flux; reflected channel modes"
    return u[:,1:half,:],v[:,1:half+1,:],w[:,1:half,:],info
end

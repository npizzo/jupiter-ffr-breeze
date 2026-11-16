module Instability2D
using Oceananigans, Statistics, LinearAlgebra, Random
using Oceananigans: TendencyCallsite
using Oceananigans.Solvers: FFTW
using Oceananigans.BoundaryConditions: fill_halo_regions!

export BandOperator, vorticity, velocity_from_curl, band_curl, spectral_curl,
       build_model, attach!, seed!, state, reflect_state, direct_channel_band

"""CPU prototype. Channel uses Fourier-x/sine-y modes, exactly as the supplied Fortran.
Periodic geometry uses the same k^2 multiplier on the full Fourier annulus.
The multiplier uses continuous physical k, while inversion uses the native discrete Laplacian.
"""
struct BandOperator
    nx::Int; ny::Int; lx::Float64; ly::Float64; channel::Bool
    k2::Matrix{Float64}; lambda::Matrix{Float64}; mask::BitMatrix
end
function BandOperator(nx,ny; lx=2pi,ly=2pi,channel=false,kmin=18.,kmax=20.,allow_zonal=false)
    nx>=8 && ny>=8 || error("Grid too small")
    0<kmin<=kmax || error("Invalid forcing band")
    dx=lx/nx;dy=ly/ny;my=channel ? 2ny : ny
    kx=[2pi/lx*(i<=nx÷2 ? i : i-nx) for i=0:nx-1]
    ky=[2pi/(channel ? 2ly : ly)*(j<=my÷2 ? j : j-my) for j=0:my-1]
    k2=[x*x+y*y for x in kx,y in ky]
    lambda=[4sin(x*dx/2)^2/dx^2+4sin(y*dy/2)^2/dy^2 for x in kx,y in ky]
    mask=BitMatrix((k2 .>= kmin^2-1e-12).&(k2 .<= kmax^2+1e-12))
    !allow_zonal && (mask[1,:].=false)
    channel && (mask[:,1].=false)
    kmax < min(pi/dx,pi/dy) || error("Band touches Nyquist")
    any(mask) || error("Empty band")
    BandOperator(nx,ny,lx,ly,channel,k2,lambda,mask)
end
function extension(o,q)
    o.channel || return copy(q)
    size(q)==(o.nx,o.ny+1)||error("Channel curl must include both wall nodes")
    a=zeros(o.nx,2o.ny);a[:,2:o.ny].=q[:,2:o.ny]
    a[:,o.ny+2:end].=-q[:,o.ny:-1:2]
    a
end
function spectral_curl(o,q,multiplier)
    a=real.(FFTW.ifft(FFTW.fft(extension(o,q)).*multiplier))
    o.channel ? hcat(zeros(o.nx),a[:,2:o.ny],zeros(o.nx)) : a
end
band_curl(o,q,a)=spectral_curl(o,q,a.*o.k2.*o.mask)
function velocity_from_curl(o,q)
    invlap=map(x->x==0 ? 0.0 : -1/x,o.lambda)
    psi=spectral_curl(o,q,invlap);dx=o.lx/o.nx;dy=o.ly/o.ny
    u=o.channel ? -(psi[:,2:end].-psi[:,1:end-1])./dy : -(circshift(psi,(0,-1)).-psi)./dy
    v=(circshift(psi,(-1,0)).-psi)./dx
    u,v
end
function vorticity(o,u,v)
    dx=o.lx/o.nx;dy=o.ly/o.ny
    if o.channel
        q=zeros(o.nx,o.ny+1)
        q[:,2:o.ny].=((v.-circshift(v,(1,0)))./dx)[:,2:o.ny].-(u[:,2:end].-u[:,1:end-1])./dy
        return q
    end
    (v.-circshift(v,(1,0)))./dx.-(u.-circshift(u,(0,1)))./dy
end
state(m)=(copy(dropdims(Array(interior(m.velocities.u));dims=3)),copy(dropdims(Array(interior(m.velocities.v));dims=3)))
function seed!(m,u,v)
    interior(m.velocities.u)[:,:,1].=u;interior(m.velocities.v)[:,:,1].=v
    fill_halo_regions!(m.velocities.u);fill_halo_regions!(m.velocities.v)
    nothing
end
function reflect_state(u,v)
    # x-face coordinates: x -> -x; x-center coordinates: i -> nx+1-i.
    nx=size(u,1);faceidx=vcat(1,collect(nx:-1:2))
    -u[faceidx,:],v[nx:-1:1,:]
end
function build_model(o; f=0.0,nu=0.0,advection=true)
    grid=RectilinearGrid(CPU(),Float64;size=(o.nx,o.ny),x=(0,o.lx),y=(0,o.ly),
                        topology=(Periodic,o.channel ? Bounded : Periodic,Flat),halo=(3,3))
    NonhydrostaticModel(grid;buoyancy=nothing,coriolis=FPlane(f=f),
                        advection=advection ? Centered(order=2) : nothing,
                        closure=nu>0 ? ScalarDiffusivity(ν=nu) : nothing,
                        timestepper=:RungeKutta3)
end
mutable struct Feedback
    op::BandOperator; amplitude::Float64; viscosity_factor::Float64; cubic_drag::Float64
    calls::Int; last_force::Any
end
function (b::Feedback)(m)
    o=b.op;u,v=state(m);q=vorticity(o,u,v)
    nu=b.viscosity_factor*sqrt(sum(abs2,q)/(o.nx*o.ny))*min(o.lx/o.nx,o.ly/o.ny)^2
    curlforce=band_curl(o,q,b.amplitude)
    nu>0 && (curlforce .+= spectral_curl(o,q,-nu.*o.lambda))
    fu,fv=velocity_from_curl(o,curlforce)
    if b.cubic_drag>0
        # Same-location interpolation on the C grid; each drag component does nonpositive work.
        vc=(v[:,1:o.ny].+(o.channel ? v[:,2:o.ny+1] : circshift(v,(0,-1))))./2
        vu=(vc.+circshift(vc,(1,0)))./2
        uc=(u.+circshift(u,(-1,0)))./2
        uv=o.channel ? hcat(uc[:,1],(uc[:,1:end-1].+uc[:,2:end])./2,uc[:,end]) : (uc.+circshift(uc,(0,1)))./2
        fu .-= b.cubic_drag .* (u .^ 2 .+ vu .^ 2) .* u
        fv .-= b.cubic_drag .* (v .^ 2 .+ uv .^ 2) .* v
    end
    interior(m.timestepper.Gⁿ.u)[:,:,1].+=fu
    interior(m.timestepper.Gⁿ.v)[:,:,1].+=fv
    b.calls+=1;b.last_force=(fu,fv)
    nothing
end
function attach!(sim,o;amplitude=1.5e-3,viscosity_factor=0.0,cubic_drag=0.0)
    viscosity_factor>=0 && cubic_drag>=0 || error("Damping must be nonnegative")
    b=Feedback(o,amplitude,viscosity_factor,cubic_drag,0,nothing)
    add_callback!(sim,b,IterationInterval(1);name=:instability_feedback,callsite=TendencyCallsite())
    b
end
"""Independent direct projection, matching add_instability without its Euler dt kick.
Fortran assumes Lx=2pi, Ly=pi. Deliberately slow: validation only.
"""
function direct_channel_band(o,q,a,kmin,kmax;allow_zonal=false)
    o.channel && o.lx≈2pi && o.ly≈pi || error("Reference geometry mismatch")
    out=zeros(size(q));nx=o.nx;ny=o.ny
    for l=1:ceil(Int,kmax),k=(allow_zonal ? 0 : 1):ceil(Int,kmax)
        kmin^2<=k*k+l*l<=kmax^2 || continue
        c=[cos(2pi*k*i/nx)*sin(pi*l*j/ny) for i=0:nx-1,j=0:ny]
        s=[sin(2pi*k*i/nx)*sin(pi*l*j/ny) for i=0:nx-1,j=0:ny]
        norm=(k==0 ? 2 : 4)/(nx*ny)
        out .+= a*(k*k+l*l)*norm.*(sum(q.*c).*c.+sum(q.*s).*s)
    end
    out[:,[1,end]].=0;out
end
end

# Additional, circulation-preserving diffusion of anticyclonic vorticity.
# No instability forcing, drag, or transform is used by this callback.
include("../viscous2d_fast512_v1/model.jl")
struct AnticycloneTerms
    op::BandOperator
    buffer::Matrix{ComplexF64}
    sign_f::Float64
    nu_A::Float64
    dt::Float64
end
@noinline function select_anticyclone!(q,s,nu_A)
    @inbounds for k in eachindex(q)
        q[k]=complex(nu_A*s*min(s*real(q[k]),0.),0.)
    end
    nothing
end
function (b::AnticycloneTerms)(m)
    o=b.op;nx,ny=o.nx,o.ny;idx=nx/o.lx;idy=ny/o.ly
    u=interior(m.velocities.u);v=interior(m.velocities.v)
    guard_cfl(u,v,nx,ny,idx,idy,b.dt)
    curl_to_buffer!(b.buffer,u,v,nx,ny,idx,idy)
    select_anticyclone!(b.buffer,b.sign_f,b.nu_A)
    add_streamfunction_tendency!(interior(m.timestepper.Gⁿ.u),interior(m.timestepper.Gⁿ.v),b.buffer,nx,ny,idx,idy)
    nothing
end
function setup_anticyclone(u,v;t=0.,f=1.,nu_A=.0004,dt=.0025)
    @assert f!=0 && nu_A>=0
    o=BandOperator(size(u,1),size(u,2);kmin=12.,kmax=16.)
    # Conservative explicit-diffusion bound, in addition to stagewise advective CFL.
    @assert dt*(.0001+nu_A)*((o.nx/o.lx)^2+(o.ny/o.ly)^2)<.25
    m=build_model(o;f,nu=.0001);seed!(m,u,v)
    m.clock.time=t;m.clock.iteration=0
    m,AnticycloneTerms(o,zeros(ComplexF64,o.nx,o.ny),sign(f),nu_A,dt)
end
function source_hashes()
    paths=["model.jl","run.jl","../viscous2d_fast512_v1/model.jl","../instability2d_v1/core.jl"]
    Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in paths)
end

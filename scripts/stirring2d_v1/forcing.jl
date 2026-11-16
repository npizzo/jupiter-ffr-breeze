include("../instability2d_long_v1/audit.jl")
include("../instability2d_long_v1/spectra.jl")

"""Independent additive pulses: same spatial spectrum, locations and expected
energy injection in both cases. Sign bias is in non-Gaussian core polarity,
not in net circulation. No amplitude feedback on the evolving velocity.
"""
mutable struct Stirring
    case::String
    f::Float64
    epsilon::Float64
    interval::Float64
    amplitude::Float64
    q::Matrix{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
    positions::MersenneTwister
    signs::MersenneTwister
    events::Int
    last_iteration::Int
    energy_quadratic::Float64
    energy_cross::Float64
    positive_events::Int
end
function Stirring(o;case,epsilon,interval=0.2,f=1.,position_seed=23017,sign_seed=23018)
    case in ("random_sign","cyclonic_sign")||error("Unknown stirring case")
    epsilon>0 && interval>0 && f!=0 || error("Invalid stirring parameters")
    # Positive coefficients give a positive centered core; zero mode absent.
    q=real.(FFTW.ifft(Float64.(o.mask)))
    u,v=velocity_from_curl(o,q);rms=sqrt(mean(u.^2 .+ v.^2))
    q./=rms;u./=rms;v./=rms
    a=sqrt(2epsilon*interval) # template kinetic energy = 1/2
    Stirring(case,f,epsilon,interval,a,q,u,v,MersenneTwister(position_seed),MersenneTwister(sign_seed),0,-1,0.,0.,0)
end
function draw_event!(s)
    ij=(rand(s.positions,0:size(s.q,1)-1),rand(s.positions,0:size(s.q,2)-1))
    coin=rand(s.signs,Bool) ? 1 : -1 # draw in BOTH cases; streams stay paired
    sg=s.case=="random_sign" ? coin : Int(sign(s.f))
    s.events+=1;s.positive_events+=(sg>0)
    ij,sg
end
function kick!(s,b,m)
    finish!(b,m) # close the RK update BEFORE adding the independent jump
    u,v=state(m);q=vorticity(b.op,u,v)
    ij,sg=draw_event!(s)
    du=(sg*s.amplitude).*circshift(s.u,ij);dv=(sg*s.amplitude).*circshift(s.v,ij)
    dq=(sg*s.amplitude).*circshift(s.q,ij)
    cross=mean(u.*du.+v.*dv);quadratic=mean(du.^2 .+dv.^2)/2
    un=u.+du;vn=v.+dv;qn=vorticity(b.op,un,vn)
    qerr=maximum(abs,qn.-q.-dq)/max(maximum(abs,qn),1e-20)
    qerr<1e-12 || error("Pulse curl mismatch")
    de=mean(un.^2 .+vn.^2)/2-mean(u.^2 .+v.^2)/2
    abs(de-cross-quadratic)<1e-12*max(mean(un.^2 .+vn.^2),1.)||error("Pulse energy closure failed")
    abs(quadratic-s.epsilon*s.interval)<1e-12||error("Expected injection normalization failed")
    dm3=mean(qn.^3)-mean(q.^3)
    b.energy[2]+=de;b.absolute_energy[2]+=abs(de)
    b.third[2]+=dm3;b.absolute_third[2]+=abs(dm3)
    s.energy_cross+=cross;s.energy_quadratic+=quadratic;s.last_iteration=m.clock.iteration
    # Exact per-event spectral ENERGY INCREMENT, not an instantaneous power.
    qh=FFTW.fft(q);dh=FFTW.fft(dq);invlam=map(x->x==0 ? 0. : 1/x,b.op.lambda)
    emode=(real.(conj.(qh).*dh).+abs2.(dh)./2).*invlam./length(q)^2
    sh=ceil.(Int,sqrt.(b.op.k2).-1e-12);shellenergy=zeros(maximum(sh)+1)
    for I in eachindex(emode);shellenergy[sh[I]+1]+=emode[I];end
    abs(sum(shellenergy)-de)<1e-11||error("Pulse spectral energy closure failed")
    maximum(abs,emode[.!b.op.mask])<1e-12||error("Out-of-band pulse work")
    seed!(m,un,vn)
    cfl=b.dt*(maximum(abs,un)/(b.op.lx/b.op.nx)+maximum(abs,vn)/(b.op.ly/b.op.ny))
    b.max_cfl=max(b.max_cfl,cfl);cfl<0.5||error("Post-pulse CFL guard exceeded: $cfl")
    # NonhydrostaticModel caches next-step stage-1 tendencies at the END of
    # the previous step. Refresh native tendencies after the jump, then let
    # Audit replace that stage's pre-jump state and add damping once.
    m.clock.stage==1||error("Pulse must be applied between full timesteps")
    Oceananigans.TimeSteppers.update_state!(m)
    b(m)
    (;ij,sg,cross,quadratic,de,dm3,shellenergy,qerr,cfl)
end
function reference_injection(o,u,v)
    q=vorticity(o,u,v);psi=spectral_curl(o,q,map(x->x==0 ? 0. : -1/x,o.lambda))
    -mean(psi.*band_curl(o,q,0.009375))
end

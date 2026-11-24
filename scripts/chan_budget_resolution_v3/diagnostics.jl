function energy(m,p)
    a=arrays(m);dv=p.dx*p.dy*p.dz;z=reshape([zcell(k,p) for k=1:p.nz],1,1,:)
    rx=.5 .* (a.rho.+circshift(a.rho,(1,0,0)))
    ry=.5 .* (a.rho.+circshift(a.rho,(0,1,0)))
    rz=.5 .* (a.rho[:,:,1:end-1].+a.rho[:,:,2:end])
    ke=.5dv*(sum(rx.*a.uf.^2)+sum(ry.*a.vf.^2)+sum(rz.*a.wf[:,:,2:end-1].^2))
    ie=dv*sum(a.rho.*CV.*a.T);pe=dv*sum(a.rho.*GRAV.*z)
    (;E=ke+ie+pe,ke,ie,pe,a)
end
mutable struct Budget
    state::NamedTuple
    q
    rates
    top
    p
end
function power(b,m)
    r=interior(density_field(m));T=interior(m.temperature)
    @. b.q=r*(T-1)*b.rates
    cooling=sum(b.q)*b.p.dx*b.p.dy*b.p.dz
    b.top .= view(T,:,:,b.p.nz)
    top=2b.p.kc*(sum(b.top)/(b.p.nx*b.p.ny)-1)/b.p.dz
    b.p.heat*(16*(b.p.flux-top)-cooling),cooling,top
end
function Budget(m,p)
    q=similar(interior(m.temperature));rates=similar(q,1,1,p.nz)
    copyto!(rates,reshape([cool_rate(k,p) for k=1:p.nz],1,1,:))
    b=Budget((;),q,rates,similar(q,p.nx,p.ny),p);P,_,_=power(b,m);e=energy(m,p);t=Float64(m.clock.time)
    b.state=(t=t,start=t,E0=e.E,power=P,heat=0.,sparse_t=t,sparse_p=P,sparse_heat=0.,mass0=mean(e.a.rho))
    b
end
function (b::Budget)(sim)
    s=b.state;t=Float64(sim.model.clock.time);P,_,_=power(b,sim.model)
    # No reference reset on checkpoint restart; all budget scalars serialize.
    b.state=merge(s,(;t,power=P,heat=s.heat+.5*(t-s.t)*(P+s.power)))
end
Oceananigans.prognostic_state(b::Budget)=b.state
Oceananigans.restore_prognostic_state!(b::Budget,s::NamedTuple)=(b.state=s;nothing)
function write_frame(m,p,out,t)
    a=arrays(m);z=[zcell(k,p) for k=1:p.nz];ks=[argmin(abs.(z.-h)) for h in (.2,.5,.9)]
    zf=(a.vf.-circshift(a.vf,(1,0,0)))/p.dx-(a.uf.-circshift(a.uf,(0,1,0)))/p.dy
    zeta=.25 .* (zf.+circshift(zf,(-1,0,0)).+circshift(zf,(0,-1,0)).+circshift(zf,(-1,-1,0)))
    jldsave(joinpath(out,@sprintf("frame_%09d.jld2",round(Int,1000t)));t,video_heights=z[ks],f=p.f,lx=p.lx,ly=p.ly,
        vorticity_slices=Float32.(zeta[:,:,ks]),tracer_slices=Float32.((Array(interior(m.tracers.dye))./a.rho)[:,:,ks]))
end

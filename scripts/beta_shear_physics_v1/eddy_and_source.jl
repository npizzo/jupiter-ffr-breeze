# Identical mixed-sign velocity seed in both cases. Balance its pressure/density
# perturbation by a reference-column displacement, then impose native vertical balance.
# This is a fresh initial-value comparison, not a claimed spun-up thermal parent.
function add_balanced_eddies!(m,p,cfg)
    a=arrays(m);psi,u,v=seed_velocity(p,cfg)
    pc=.25 .* (psi[:,1:end-1].+psi[:,2:end].+circshift(psi[:,1:end-1],(-1,0)).+circshift(psi[:,2:end],(-1,0)))
    ue=.5 .* (u.+circshift(u,(-1,0)));ve=.5 .* (v[:,1:end-1].+v[:,2:end])
    Q=p.f.*pc.-.5 .* (ue.^2 .+ve.^2 .+(2pi/4)^2*5 .*pc.^2);Q.-=mean(Q)
    z=[zcell(k,p) for k=1:p.nz]
    P=Array(interior(m.dynamics.pressure));rho=copy(a.rho)
    for k=1:p.nz,j=1:p.ny,i=1:p.nx
        ref=reference(z[k]);shift=reference(z[k]-Q[i,j]/GRAV)
        P[i,j,k]*=shift.p/ref.p;rho[i,j,k]*=shift.rho/ref.rho
    end
    # Candidate tilted-rotation balance correction: replace the jet-only
    # effective vertical gravity g-h*U by that of the combined jet+eddy flow.
    # Preserve prescribed face velocities; normalize total mass afterwards.
    if get(cfg,"eddy_tilt_correction",false)
        ue3=reshape(ue,p.nx,p.ny,1)
        rho .*= (GRAV .-p.fh.*a.u)./(GRAV .-p.fh.*(a.u.+ue3))
    end
    scale=mean(a.rho)/mean(rho);rho.*=scale;P.*=scale
    U=a.uf.+reshape(u,p.nx,p.ny,1);V=a.vf.+reshape(v,p.nx,p.ny+1,1)
    uc=.5 .* (U.+circshift(U,(-1,0,0)))
    Pr=vec(Array(interior(m.dynamics.reference_state.pressure)));rr=vec(Array(interior(m.dynamics.reference_state.density)))
    for k=p.nz-1:-1:1
        rz=.5 .* (rho[:,:,k].+rho[:,:,k+1]);rref=.5*(rr[k]+rr[k+1])
        mz=.5 .* (rho[:,:,k].*uc[:,:,k].+rho[:,:,k+1].*uc[:,:,k+1])
        P[:,:,k].=P[:,:,k+1].-(Pr[k+1]-Pr[k].+p.dz.*(-GRAV.*(rz.-rref).+p.fh.*mz))
    end
    @assert minimum(P)>0 && minimum(rho)>0
    set!(m;ρ=rho,θ=P.^.6 ./rho,u=U,v=V,w=0.,qᵗ=0.)
    Oceananigans.TimeSteppers.update_state!(m)
    @assert abs(mean(Array(interior(m.dynamics.dry_density)))/mean(a.rho)-1)<1e-12
end

# Explicit source depletion bound: at most 10% of rho*theta per outer step
# according to the current complete heating source. This is a safety bound,
# not a claim of unconditional RK positivity; stage guards remain active.
@kernel function source_fraction_rates!(out,g,F,p)
    i,j,k=@index(Global,NTuple)
    q=heat_source(i,j,k,g,nothing,F,p)
    out[i,j,k]=max(0.,-q)/max(2.5F.ρᵈ[i,j,k]*F.T[i,j,k],1e-30)
end

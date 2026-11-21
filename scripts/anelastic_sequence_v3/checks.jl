include("model.jl")
root=ARGS[1];mkpath(root)
c=TOML.parsefile("configs/anelastic-sequence-v3/baseline.toml")
small=merge(c,Dict("nx"=>8,"ny"=>8,"nz"=>32))
results=Dict{String,Any}()
function report!(name,d)
    results[name]=d;record(joinpath(root,"CHECKS.toml"),results)
    println("CHECK ",name," ",d);flush(stdout)
end
for (label,caplog) in (("weak",log(1.01)),("strong",.25))
    cc=merge(small,Dict("cap_log_theta"=>caplog));m,p=build(cc;thermal=false,rest=true,inviscid=true)
    start=fields(m).theta;evolve!(m,.2)
    inv=invariants(m,p);err=maximum(abs,fields(m).theta.-start)
    @assert inv.speed<1e-10 && inv.div<1e-10 && err<1e-10
    report!(label*"_rest",Dict("status"=>"PASS","speed"=>inv.speed,"theta_error"=>err,"mass_divergence"=>inv.div))
end
# Closed integral heat test at actual weak-cap stratification; integrate the known cooling source.
m,p=build(small);inv0=invariants(m,p);prev=Ref((t=0.,source=0.));acc=Ref(0.)
function source(m,p)
    a=fields(m);z=[(k-.5)*p.dz for k=1:p.nz];target=reshape(theta_reference.(z,Ref(p)),1,1,:);rate=reshape(shape.(z,Ref(p)),1,1,:)/p.tau
    p.length^2*p.flux-sum(a.rho.*(a.theta.-target).*rate)*(p.length/p.nx)*(p.length/p.ny)*p.dz
end
prev[]=(t=0.,source=source(m,p))
function heat_cb(sim)
    t=Float64(sim.model.clock.time);s=source(sim.model,p);old=prev[]
    acc[]+=(t-old.t)*(s+old.source)/2;prev[]=(t=t,source=s)
end
evolve!(m,.2;callback=heat_cb)
err=abs(invariants(m,p).mass-inv0.mass-acc[])/(p.length^2*p.flux*.2)
@assert err<1e-5
report!("thermal_integral",Dict("status"=>"PASS","error_over_imposed_input"=>err))
# Proper reflection x -> L-x, including staggered x-face indices and f -> -f.
mp,pp=build(small);mm,pm=build(merge(small,Dict("f"=>-small["f"])))
set!(mm;θ=reverse(fields(mp).theta;dims=1),qᵗ=0.,u=0.,v=0.,w=0.)
evolve!(mp,.4);evolve!(mm,.4)
a=fields(mp);b=fields(mm);perm=vcat(1,collect(pp.nx:-1:2))
errs=Dict("u"=>maximum(abs,b.u.+a.u[perm,:,:]),"v"=>maximum(abs,b.v.-reverse(a.v;dims=1)),"w"=>maximum(abs,b.w.-reverse(a.w;dims=1)),"theta"=>maximum(abs,b.theta.-reverse(a.theta;dims=1)))
@assert maximum(values(errs))<1e-10
report!("reflected_f",Dict("status"=>"PASS","absolute_errors"=>errs,"scope"=>"Full heated weak-cap model; x-reflected temperature and staggered velocities, opposite vertical f"))
# Linear internal-wave eigenmode for the actual strong-cap hydrostatic reference.
# For w=psi/sqrt(rho), A=-dzz+s²/4-s'/2 and omega²(A+kh²)=f² A+kh² N².
# This independent continuum eigenproblem is sampled on a vertical finite-difference mesh.
wc=merge(c,Dict("nx"=>32,"ny"=>8,"nz"=>64,"length"=>.5,"cap_log_theta"=>.25))
m,p=build(wc;thermal=false,rest=true,inviscid=true)
a=fields(m);rho=vec(a.rho);Tref=vec(Array(interior(p.Tref)));nz=p.nz;dz=p.dz;dx=p.length/p.nx;kh=2pi/p.length
rf=vcat(rho[1],0.5 .*(rho[1:end-1].+rho[2:end]),rho[end]);Tf=0.5 .*(Tref[1:end-1].+Tref[2:end]);zz=collect(1:nz-1).*dz
s=zeros(nz-1);sp=similar(s);N2=similar(s)
for k=1:nz-1
    x=(zz[k]-p.dc)/(p.height-p.dc)
    q=0<x<1 ? p.caplog*6x*(1-x)/(p.height-p.dc) : 0.
    qp=0<x<1 ? p.caplog*(6-12x)/(p.height-p.dc)^2 : 0.
    dT=Tf[k]*q-p.g/p.cp
    s[k]=-p.g*(1/p.R-1/p.cp)/Tf[k]-q
    sp[k]=p.g*(1/p.R-1/p.cp)*dT/Tf[k]^2-qp
    N2[k]=p.g*q
end
A=Matrix(SymTridiagonal(2/dz^2 .+ s.^2/4 .-sp/2,fill(-1/dz^2,nz-2)))
B=A+kh^2*I;C=p.f^2*A+kh^2*Diagonal(N2)
ev=eigen(Symmetric(C),Symmetric(B));omega=sqrt(ev.values[end]);psi=ev.vectors[:,end]
mode=vcat(0.,psi./sqrt.(rf[2:end-1]),0.);mode./=maximum(abs,mode)
amp=1e-6;kd=2sin(kh*dx/2)/dx
mz=[amp*rf[k]*mode[k]*cos(kh*(i-.5)*dx) for i=1:p.nx,j=1:p.ny,k=1:nz+1]
mx=[-amp*(rf[k+1]*mode[k+1]-rf[k]*mode[k])/dz/kd*sin(kh*(i-1)*dx) for i=1:p.nx,j=1:p.ny,k=1:nz]
set!(m;ρu=mx,ρv=0.,ρw=mz)
shape3=[mode[k]*cos(kh*(i-.5)*dx) for i=1:p.nx,j=1:p.ny,k=1:nz+1];normmode=sum(shape3.^2)
ts=Float64[];ys=Float64[]
function wave_cb(sim)
    push!(ts,sim.model.clock.time);push!(ys,sum(fields(sim.model).w.*shape3)/normmode/amp)
end
evolve!(m,1.5*2pi/omega;dt=.02,callback=wave_cb,interval=5)
# The end time need not align with the recording cadence; fit all sampled points.
trials=collect(range(.75omega,1.25omega;length=401));errors=Float64[]
for om in trials
    X=hcat(cos.(om.*ts),sin.(om.*ts),ones(length(ts)));coeff=X\ys;push!(errors,sqrt(mean((X*coeff.-ys).^2)))
end
best=argmin(errors);measured=trials[best];rel=abs(measured/omega-1)
record(joinpath(root,"WAVE.toml"),Dict("predicted_frequency"=>omega,"measured_frequency"=>measured,"relative_error"=>rel,"fit_rms"=>errors[best],"sample_count"=>length(ts)))
@assert rel<.05 && errors[best]<.15
report!("strong_cap_wave",Dict("status"=>"PASS","predicted_frequency"=>omega,"measured_frequency"=>measured,"relative_error"=>rel,"fit_rms"=>errors[best]))
record(joinpath(root,"PASS.toml"),Dict("status"=>"PASS","model_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"model.jl")))),"scope"=>"CPU weak/strong reference rest, weak thermal integral, reflected f, strong-cap wave. Not GPU/restart/resolution/steady-convection certification.","automatic_production_release"=>false))
println("ANELASTIC_SEQUENCE_CPU_PASS");flush(stdout)

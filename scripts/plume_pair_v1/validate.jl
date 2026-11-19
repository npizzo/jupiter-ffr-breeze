include("run.jl")
root=ARGS[1];mkpath(root)
cfg=TOML.parsefile("configs/plume-pair-v1/stable.toml")
c=merge(cfg,Dict("nx"=>16,"ny"=>16,"nz"=>40,"dt"=>.5,"stop_seconds"=>2.,"output_seconds"=>1.,"wall_time_limit"=>600.))
# Discrete centered curl of constant flow and analytic sinusoid.
n=16;h=cfg["lx"]/n;x=((1:n).-0.5).*h
u=zeros(n,n);v=repeat(reshape(sin.(2pi.*x./cfg["lx"]),:,1),1,n)
expected=repeat(reshape(cos.(2pi.*x./cfg["lx"]).*sin(2pi/n)./h,:,1),1,n)
@assert maximum(abs,pcurl(u,v,h).-expected)<1e-18
@assert maximum(abs,pcurl(ones(n,n),ones(n,n),h))==0
rest=make_plume(merge(c,Dict("pulse_kelvin"=>0.)),CPU())
@assert all(iszero,Array(interior(rest.model.velocities.w)))
@assert rest.p70>0 && all(>(0),rest.rho)
println("PLUME_ANALYTIC_AND_REFERENCE_PASS");flush(stdout)
a=run_plume(c,joinpath(root,"first");arch=CPU())
cp=only(filter(f->endswith(f,".jld2"),readdir(joinpath(root,"first/checkpoints");join=true)))
b=run_plume(merge(c,Dict("stop_seconds"=>4.)),joinpath(root,"resumed");arch=CPU(),pickup=cp)
d=run_plume(merge(c,Dict("stop_seconds"=>4.)),joinpath(root,"continuous");arch=CPU())
for key in (:ρu,:ρv,:ρw)
    aa=Array(interior(getproperty(b.q.model.momentum,key)));bb=Array(interior(getproperty(d.q.model.momentum,key)))
    @assert maximum(abs,aa.-bb)<=1e-11+1e-8maximum(abs,bb)
end
@assert b.b.checks>0 && d.b.checks>=8
neutral=run_plume(merge(c,Dict("cap_n2"=>0.)),joinpath(root,"neutral");arch=CPU())
@assert neutral.b.checks>=4
open(joinpath(root,"PASS.toml"),"w") do io
    TOML.print(io,Dict("status"=>"PASS","checks"=>"analytic curl; hydrostatic reference; native layer RHS partition and SSPRK closure in both controls; serialized velocity restart","source_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"run.jl")))),"budget_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"budget.jl")))),"scope"=>"CPU implementation check; GPU budget checks remain active throughout production; not resolution convergence"))
end
println("PLUME_PAIR_IMPLEMENTATION_PASS");flush(stdout)

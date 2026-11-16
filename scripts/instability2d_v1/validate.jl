include("core.jl")
using .Instability2D, Oceananigans, Statistics, LinearAlgebra, Random, TOML, Test
out=ARGS[1];ispath(out)&&error("Refusing to overwrite validation output");mkpath(out)
metrics=Dict{String,Any}();rng=MersenneTwister(1701)
relative(a,b)=norm(a-b)/max(norm(b),1e-30)
@testset "Fortran projection, native curl and symmetry" begin
    for zonal in (false,true)
        o=BandOperator(64,32;channel=true,ly=pi,kmin=3.,kmax=6.,allow_zonal=zonal)
        q=randn(rng,64,33);q[:,[1,end]].=0
        a=band_curl(o,q,.0015);ref=direct_channel_band(o,q,.0015,3,6;allow_zonal=zonal)
        err=relative(a,ref);metrics["fortran_projection_zonal_$(zonal)"]=err;@test err<1e-12
        fu,fv=velocity_from_curl(o,a)
        @test relative(vorticity(o,fu,fv),a)<1e-12
        div=(circshift(fu,(-1,0)).-fu)./(o.lx/o.nx).+(fv[:,2:end].-fv[:,1:end-1])./(o.ly/o.ny)
        @test maximum(abs,div)<1e-12
    end
end
function advance(o,u,v;f=0.,a=.0015,dt=.001,t=.1,advection=true)
    m=build_model(o;f,advection);seed!(m,u,v)
    sim=Simulation(m;Δt=dt,stop_time=t,verbose=false)
    b=attach!(sim,o;amplitude=a);run!(sim)
    @test b.calls>=3m.clock.iteration
    state(m),b.calls
end
@testset "Oceananigans integration" begin
    for channel in (false,true)
        o=BandOperator(32,32;channel,ly=channel ? pi : 2pi,kmin=3.,kmax=6.)
        y=channel ? collect(0:32).*pi/32 : collect(0:31).*2pi/32
        q=[cos(3*2pi*i/32)*sin(4*yy) for i=0:31,yy=y]
        u,v=velocity_from_curl(o,q)
        (uf,vf),calls=advance(o,u,v;advection=false)
        expected=exp(.0015*25*.1)
        err=relative(vorticity(o,uf,vf),expected.*q)
        metrics["mode_growth_error_channel_$(channel)"]=err;@test err<1e-10
        @test abs(sum(u.*u)+sum(v.*v))>0
        # In-band velocity acceleration is the same scalar growth as vorticity.
        fu,fv=velocity_from_curl(o,band_curl(o,q,.0015))
        @test sum(u.*fu)+sum(v.*fv)>0
        qo=[cos(2pi*i/32)*sin(yy) for i=0:31,yy=y]
        @test norm(band_curl(o,qo,.0015))<1e-12
        # A nonlinear multi-mode field, with zero mean velocity.
        q .+= 0.3 .* [sin(2*2pi*i/32)*sin(3*yy)+0.2*cos(5*2pi*i/32)*sin(yy) for i=0:31,yy=y]
        u,v=velocity_from_curl(o,q);u.*=.01;v.*=.01
        (up,vp),_=advance(o,u,v;f=1.)
        (um,vm),_=advance(o,u,v;f=-1.)
        (uz,vz),_=advance(o,u,v;f=0.)
        e=max(relative(up,uz),relative(vp,vz),relative(um,uz),relative(vm,vz))
        metrics["f_independence_channel_$(channel)"]=e;@test e<1e-8
        ur,vr=reflect_state(u,v);(ua,va),_=advance(o,ur,vr;f=0.)
        ue,ve=reflect_state(uz,vz)
        er=max(relative(ua,ue),relative(va,ve));metrics["reflection_error_channel_$(channel)"]=er;@test er<1e-10
    end
end
metrics["status"]="PASS";metrics["scope"]="Small CPU operator/integration checks only; no claim about turbulent sign selection"
metrics["compute_location"]="Unity";metrics["oceananigans_version"]=string(pkgversion(Oceananigans))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
println("INSTABILITY2D_VALIDATION_PASS");TOML.print(stdout,metrics;sorted=true)

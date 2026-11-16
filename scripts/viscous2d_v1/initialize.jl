include("../stirring2d_v1/forcing.jl")
using JLD2,SHA,Test
out=ARGS[1];ispath(out)&&error("Refusing existing source");mkpath(out)
o=BandOperator(256,256;kmin=12.,kmax=16.,allow_zonal=false)
rng=MersenneTwister(62041)
q=spectral_curl(o,randn(rng,256,256),(o.k2.>=12^2).&(o.k2.<=16^2))
u,v=velocity_from_curl(o,q);factor=.2/sqrt(mean(u.^2 .+v.^2));u.*=factor;v.*=factor
q=vorticity(o,u,v);E=mean(u.^2 .+v.^2)/2
s=Stirring(o;case="random_sign",epsilon=reference_injection(o,u,v)*(.0002/.009375))
net=0.0002 .* o.k2 .- 0.0001 .* o.lambda
@testset "Viscosity-only setup" begin
 @test isapprox(sqrt(2E),.2;rtol=1e-12)
 @test abs(mean(u))+abs(mean(v))+abs(mean(q))<1e-12
 @test maximum(abs,FFTW.fft(q)[(o.k2.<12^2).|(o.k2.>16^2)])<1e-9
 @test all(net[o.mask].>0)
 @test sort(unique(mod.(81 .* collect(1:80),80)))==collect(0:79)
 @test isapprox(s.amplitude^2/2,s.epsilon*.2;rtol=1e-12)
 @test Audit(o;alpha=0.).alpha==0
end
jldsave(joinpath(out,"initial.jld2");u,v,t=0.,seed=62041,initial_band=[12,16])
qh=FFTW.fft(q);invlam=map(x->x==0 ? 0. : 1/x,o.lambda);ek=abs2.(qh).*invlam
invk=map(x->x==0 ? 0. : 1/sqrt(x),o.k2);L=2pi*sum(ek.*invk)/sum(ek)
params=Dict("initial_urms"=>sqrt(2E),"initial_energy"=>E,"initial_integral_wavelength"=>L,
 "initial_turnover_L_over_U"=>L/.2,"initial_energy_input_time"=>E/s.epsilon,
 "stochastic_expected_injection"=>s.epsilon,"pulse_velocity_ratio"=>s.amplitude/.2,
 "net_linear_growth_range"=>collect(extrema(net[o.mask])),"nu"=>.0001,"drag"=>0.,"a"=>.0002,
 "initial_sha256"=>bytes2hex(sha256(read(joinpath(out,"initial.jld2")))))
open(joinpath(out,"parameters.toml"),"w") do io;TOML.print(io,params;sorted=true);end
TOML.print(stdout,params;sorted=true)

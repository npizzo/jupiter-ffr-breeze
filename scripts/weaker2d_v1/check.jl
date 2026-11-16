include("../stirring2d_v1/forcing.jl")
using Test,JLD2,SHA
source,out=ARGS;d=load(source);u,v=d["u"],d["v"]
o=BandOperator(128,128;kmin=6.,kmax=8.,allow_zonal=false)
eps0=reference_injection(o,u,v);eps=.2eps0
s=Stirring(o;case="random_sign",epsilon=eps)
q=vorticity(o,u,v);U=sqrt(mean(u.^2 .+v.^2));E=U^2/2
net=0.001875 .* o.k2 .- 0.001 .* o.lambda
metrics=Dict{String,Any}()
@testset "Weaker forcing and sampling" begin
 @test all(net[o.mask].>0)
 @test all((0.0009375 .* o.k2 .- 0.001 .* o.lambda)[o.mask].<0)
 @test isapprox(s.amplitude/sqrt(2eps0*.2),sqrt(.2);rtol=1e-14)
 @test sort(unique(mod.(201 .* collect(1:40),40)))==collect(0:39)
 @test norm(band_curl(o,q,.001875)-.2band_curl(o,q,.009375))/norm(q)<1e-12
 @test isapprox(s.amplitude^2/2,eps*.2;rtol=1e-14)
 @test abs(mean(s.q))<1e-12
 @test norm(vorticity(o,s.u,s.v)-s.q)/norm(s.q)<1e-12
 # Net modal growth includes native-grid viscous damping. No drag change.
end
qh=FFTW.fft(q);invlam=map(x->x==0 ? 0. : 1/x,o.lambda);ek=abs2.(qh).*invlam
invk=map(x->x==0 ? 0. : 1/sqrt(x),o.k2);L=2pi*sum(ek.*invk)/sum(ek)
metrics["initial_urms"]=U;metrics["initial_integral_wavelength"]=L;metrics["initial_turnover_time_L_over_U"]=L/U
metrics["old_reference_energy_input_time"]=E/eps0;metrics["new_reference_energy_input_time"]=E/eps
metrics["old_pulse_velocity_ratio"]=sqrt(2eps0*.2)/U;metrics["new_pulse_velocity_ratio"]=s.amplitude/U
metrics["net_linear_growth_range"]=extrema(net[o.mask]) |> collect
metrics["old_expected_injection"]=eps0;metrics["new_expected_injection"]=eps
metrics["source_checkpoint_sha256"]=bytes2hex(sha256(read(source)))
open(out,"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true)

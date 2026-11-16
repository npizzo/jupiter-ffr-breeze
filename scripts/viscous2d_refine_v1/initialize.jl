include("../stirring2d_v1/forcing.jl")
using JLD2,SHA,Test
source,out=ARGS
ispath(out)&&error("Refusing existing output");mkpath(out)
d=load(source);uc,vc=d["u"],d["v"];@assert size(uc)==size(vc)==(256,256) && d["t"]==0
oc=BandOperator(256,256;kmin=12.,kmax=16.,allow_zonal=false)
of=BandOperator(512,512;kmin=12.,kmax=16.,allow_zonal=false)
qc=vorticity(oc,uc,vc);ch=FFTW.fft(qc)./256^2;fh=zeros(ComplexF64,512,512)
# Curl lives at horizontal cell corners on both grids, with the same origin.
# Preserve each signed Fourier coefficient. No random redraw or RMS rescaling.
for kx=-16:16,ky=-16:16
    if 12^2<=kx^2+ky^2<=16^2
        fh[mod(kx,512)+1,mod(ky,512)+1]=ch[mod(kx,256)+1,mod(ky,256)+1]
    end
end
qf=real.(FFTW.ifft(fh.*512^2));u,v=velocity_from_curl(of,qf)
recovered=vorticity(of,u,v);E0=mean(uc.^2 .+vc.^2)/2;E1=mean(u.^2 .+v.^2)/2
restriction_error=maximum(abs,qf[1:2:end,1:2:end].-qc)/maximum(abs,qc)
div=(circshift(u,(-1,0)).-u)./(2pi/512).+(circshift(v,(0,-1)).-v)./(2pi/512)
@testset "Preserve the original vorticity realization" begin
 @test restriction_error<1e-12
 @test norm(recovered-qf)/norm(qf)<1e-12
 @test maximum(abs,div)<1e-12
 @test abs(mean(u))+abs(mean(v))+abs(mean(qf))<1e-12
 @test abs(E1/E0-1)<.01
 @test maximum(abs,FFTW.fft(qf)[(of.k2.<12^2).|(of.k2.>16^2)])<1e-9
end
jldsave(joinpath(out,"initial.jld2");u,v,t=0.,source,source_sha256=bytes2hex(sha256(read(source))))
metadata=Dict("method"=>"Preserve normalized corner-vorticity Fourier coefficients, then invert the fine native curl; no amplitude renormalization",
 "coarse_initial_energy"=>E0,"fine_initial_energy"=>E1,"relative_initial_energy_change"=>E1/E0-1,
 "coarse_urms"=>sqrt(2E0),"fine_urms"=>sqrt(2E1),"relative_vorticity_restriction_error"=>restriction_error,
 "max_divergence"=>maximum(abs,div),"source_sha256"=>bytes2hex(sha256(read(source))),
 "initial_sha256"=>bytes2hex(sha256(read(joinpath(out,"initial.jld2")))))
open(joinpath(out,"initialization.toml"),"w") do io;TOML.print(io,metadata;sorted=true);end
TOML.print(stdout,metadata;sorted=true)

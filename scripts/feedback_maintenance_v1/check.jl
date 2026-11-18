include("model.jl")
using Test
source,out=ARGS
ispath(out)&&error("Output exists");mkpath(out)
d=load(source);check_source_checkpoint(d);u0,v0=d["u"],d["v"]
origin_hash=bytes2hex(sha256(read(source)))
rel(a,b)=sqrt(sum(abs2,a.-b)/max(sum(abs2,b),eps()))
function reference_maintenance!(m,o;remove_low)
    u,v=state(m);q=vorticity(o,u,v);s=maintenance_sources(o,q;remove_low)
    fu,fv=velocity_from_curl(o,s.linear.+s.quadratic.+s.cubic)
    interior(m.timestepper.Gⁿ.u)[:,:,1].+=fu;interior(m.timestepper.Gⁿ.v)[:,:,1].+=fv
end
function evolve_maintenance(u,v;case="full",steps=20,t=400.,dt=.001,reference=false)
    m,b=setup_maintenance(u,v;case,t,dt)
    sim=Simulation(m;Δt=dt,stop_iteration=steps,verbose=false)
    cb=reference ? (m->reference_maintenance!(m,b.base.op;remove_low=b.remove_low)) : b
    add_callback!(sim,cb,IterationInterval(1);name=:maintenance,callsite=TendencyCallsite());run!(sim)
    (;m,b,sim)
end
# Independent centered finite-volume flux reconstruction used by the old-output audit.
function reconstructed_advection_curl(o,u,v)
    dx=o.lx/o.nx;dy=o.ly/o.ny
    uc=.5 .* (u.+circshift(u,(-1,0)));vc=.5 .* (v.+circshift(v,(0,-1)))
    cross=.25 .* (u.+circshift(u,(0,1))) .* (v.+circshift(v,(1,0)))
    au=-(uc.^2 .-circshift(uc.^2,(1,0)))./dx .-(circshift(cross,(0,-1)).-cross)./dy
    av=-(circshift(cross,(-1,0)).-cross)./dx .-(vc.^2 .-circshift(vc.^2,(0,1)))./dy
    vorticity(o,au,av)
end
metrics=Dict{String,Any}()
@testset "Matched maintenance ablation" begin
    m,b=setup_maintenance(u0,v0;case="full");o=b.base.op
    native_reconstruction_error=rel(native_advection_curl(m,o),reconstructed_advection_curl(o,u0,v0))
    @test native_reconstruction_error<1e-10
    q=vorticity(o,u0,v0);sfull=maintenance_sources(o,q);scut=maintenance_sources(o,q;remove_low=true)
    @test sfull.cubic==scut.cubic && sfull.linear==scut.linear && sfull.viscous==scut.viscous
    @test !any(o.mask[b.low]) && !any(o.mask[1,:])
    diff=FFTW.fft(sfull.quadratic.-scut.quadratic)
    @test maximum(abs,diff[.!b.low])<1e-10*maximum(abs,diff)
    @test maximum(abs,FFTW.fft(scut.quadratic)[b.low])<1e-10*maximum(abs,diff)
    @test abs(mean(u0))<1e-12 && abs(mean(v0))<1e-12
    curlerrs=Float64[];dt=.001
    for case in MAINTENANCE_CASES
        mm,bb=setup_maintenance(u0,v0;case,dt)
        fill!(interior(mm.timestepper.Gⁿ.u),0.);fill!(interior(mm.timestepper.Gⁿ.v),0.);bb(mm)
        fu=copy(dropdims(Array(interior(mm.timestepper.Gⁿ.u));dims=3));fv=copy(dropdims(Array(interior(mm.timestepper.Gⁿ.v));dims=3))
        ss=maintenance_sources(o,q;remove_low=bb.remove_low);expected=ss.linear.+ss.quadratic.+ss.cubic
        err=rel(vorticity(o,fu,fv),expected);push!(curlerrs,err);@test err<1e-10
        @test abs(mean(fu))<1e-12 && abs(mean(fv))<1e-12
        @test abs(mean(expected))<1e-12
        psi=spectral_curl(o,q,map(x->x==0 ? 0. : -1/x,o.lambda))
        @test isapprox(mean(u0.*fu.+v0.*fv),-mean(psi.*expected);rtol=1e-10,atol=1e-12)
        qq,di=maintenance_diagnostics(mm,bb)
        @test di["linear_lowk_energy_rate"]≈0 atol=1e-12
        @test abs(di["advection_energy_rate"])<1e-10
        @test isapprox(di["viscous_energy_rate"],-.0001mean(q.^2);rtol=1e-10,atol=1e-12)
        bb.remove_low && (@test abs(di["quadratic_lowk_energy_rate"])<1e-12)
        @test mean(q.*ss.cubic)<0
    end
    # Independent native G capture at a real TendencyCallsite, before closure addition.
    mm,bb=setup_maintenance(u0,v0;case="full",dt);nativeerr=Ref(NaN);observed=Ref(false)
    sim=Simulation(mm;Δt=dt,stop_iteration=1,verbose=false)
    function native_capture(m)
        if !observed[]
            u,v=state(m);qq=vorticity(o,u,v)
            gu=copy(dropdims(Array(interior(m.timestepper.Gⁿ.u));dims=3));gv=copy(dropdims(Array(interior(m.timestepper.Gⁿ.v));dims=3))
            expected=native_advection_curl(m,o).+spectral_curl(o,qq,-.0001 .*o.lambda)
            nativeerr[]=rel(vorticity(o,gu,gv),expected);observed[]=true
        end
        bb(m)
    end
    add_callback!(sim,native_capture,IterationInterval(1);name=:capture,callsite=TendencyCallsite());run!(sim)
    @test observed[] && nativeerr[]<1e-9
    # Full wrapper retains the old forcing implementation exactly.
    oldm,oldb=setup_polynomial(u0,v0;case="biased_plus",t=400.,dt)
    oldsim=Simulation(oldm;Δt=dt,stop_iteration=20,verbose=false)
    add_callback!(oldsim,oldb,IterationInterval(1);name=:old,callsite=TendencyCallsite());run!(oldsim)
    full=evolve_maintenance(u0,v0;dt);uf,vf=state(full.m);uo,vo=state(oldm)
    olderr=max(rel(uf,uo),rel(vf,vo));@test olderr<1e-12
    ab=evolve_maintenance(u0,v0;case="no_low_quadratic",dt)
    ref=evolve_maintenance(u0,v0;case="no_low_quadratic",dt,reference=true)
    ua,va=state(ab.m);ur,vr=state(ref.m);referr=max(rel(ua,ur),rel(va,vr));@test referr<1e-10
    half=evolve_maintenance(u0,v0;case="no_low_quadratic",dt,steps=9);u,v=state(half.m)
    path=joinpath(out,"restart.jld2");jldsave(path;u,v,t=Float64(half.m.clock.time))
    dd=load(path);res=evolve_maintenance(dd["u"],dd["v"];case="no_low_quadratic",dt,steps=11,t=dd["t"])
    u,v=state(res.m);restart=max(rel(u,ua),rel(v,va));@test restart<1e-10
    fine=evolve_maintenance(u0,v0;case="no_low_quadratic",dt=dt/2,steps=40)
    u,v=state(fine.m);dterr=max(rel(u,ua),rel(v,va));@test dterr<1e-3
    ab.sim.stop_iteration=70;GC.gc();perf=@timed run!(ab.sim)
    @test perf.bytes/50<65536 && perf.time/50<.25
    # Exact output schedule and adaptive step at the mature state.
    pilots=Dict{String,Any}()
    for case in MAINTENANCE_CASES
        pm,pb=setup_maintenance(u0,v0;case);ps=Simulation(pm;Δt=.0025,stop_time=400.5,verbose=false)
        maintenance_step!(ps,pb);saved=Float64[]
        add_callback!(ps,pb,IterationInterval(1);name=:terms,callsite=TendencyCallsite())
        add_callback!(ps,s->maintenance_step!(s,pb),IterationInterval(10);name=:adaptive)
        add_callback!(ps,s->push!(saved,Float64(s.model.clock.time)),TimeInterval(.25);name=:samples)
        run!(ps);pq,pd=maintenance_diagnostics(pm,pb)
        @test pm.clock.time≈400.5 atol=1e-9
        @test all(x->any(abs.(saved.-x).<1e-8),[400.25,400.5])
        @test all(isfinite,pq) && abs(mean(pq))<1e-12
        pilots[case]=Dict("t"=>pm.clock.time,"dt"=>ps.Δt,"energy"=>pd["energy"],"lowk_energy"=>pd["energy_k_le_2"])
    end
    merge!(metrics,Dict("native_advection_diffusion_curl_relative_error"=>nativeerr[],"old_output_flux_reconstruction_relative_error"=>native_reconstruction_error,"force_curl_max_relative_error"=>maximum(curlerrs),
      "old_model_relative_error"=>olderr,"independent_ablation_relative_error"=>referr,"restart_relative_error"=>restart,
      "half_dt_relative_difference"=>dterr,"seconds_per_step"=>perf.time/50,"bytes_per_step"=>perf.bytes/50,"pilots"=>pilots))
end
merge!(metrics,Dict("status"=>"PASS","source_checkpoint_sha256"=>origin_hash,"code_sha256"=>maintenance_hashes(),"target_limit"=>450.))
open(joinpath(out,"PASS.toml"),"w") do io;TOML.print(io,metrics;sorted=true);end
TOML.print(stdout,metrics;sorted=true);flush(stdout)
